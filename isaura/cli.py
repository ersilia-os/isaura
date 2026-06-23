import contextlib, datetime, os, sys
import rich_click as click
import rich_click.rich_click as rc
from isaura.manage import (
  IsauraCopy,
  IsauraMolRemover,
  IsauraWriter,
  IsauraReader,
  IsauraInspect,
  IsauraPull,
  IsauraPush,
  IsauraStat,
)
from isaura.helpers import (
  DEFAULT_BUCKET_NAME,
  DEFAULT_PRIVATE_BUCKET_NAME,
  logger,
  console,
  inspect_table,
  inspect_table_cloud,
  make_table,
  show_figlet,
  run_docker_compose,
  get_engine_status,
)

click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_ARGUMENTS = True
rc.USE_RICH_MARKUP = True
rc.SHOW_ARGUMENTS = True
rc.COLOR_SYSTEM = "truecolor"
rc.STYLE_OPTION = "bold magenta"
rc.STYLE_COMMAND = "bold green"
rc.STYLE_METAVAR = "italic yellow"
rc.STYLE_SWITCH = "underline cyan"
rc.STYLE_USAGE = "bold blue"
rc.STYLE_OPTION_DEFAULT = "dim italic"


def apply_opts(*opts):
  def _wrap(f):
    for opt in reversed(opts):
      f = opt(f)
    return f

  return _wrap


def _check_model_filename(model_id, filepath):
  """Warn if a model ID found in the filename doesn't match model_id."""
  import re
  if not filepath:
    return
  basename = os.path.basename(filepath)
  found = re.findall(r'eos[a-z0-9]{4}', basename.lower())
  for match in found:
    if match != model_id.lower():
      console.print(
        f"[red]Filename and --model-id do not match:[/] "
        f"[bold]{match}[/bold] (filename) vs [bold]{model_id}[/bold] (--model-id)."
      )
      sys.exit(1)


def _resolve_version(model_id, bucket, store):
  """Return the latest version stored for a model in a bucket, or 'v1' as fallback."""
  try:
    keys = [obj["Key"] for obj in store.list_keys(bucket, f"{model_id}/")]
    versions = {k.split("/")[1] for k in keys if len(k.split("/")) >= 2 and k.split("/")[1].startswith("v")}
    if not versions:
      return "v1"
    return max(versions, key=lambda v: int(v[1:]) if v[1:].isdigit() else 0)
  except Exception:
    return "v1"


@click.group()
def cli():
  pass


if "--help" in sys.argv or "-h" in sys.argv or len(sys.argv) == 1:
  show_figlet()
opt_model = click.option("--model-id", "-m", "model", required=True, help="Ersilia model id (eosxxxx)")
opt_model_opt = click.option("--model-id", "-m", "model", required=False, default=None, help="Ersilia model id (eosxxxx)")
opt_version = click.option("--version", "-v", default="v1", show_default=True, help="Model version")
opt_project = click.option(
  "--project-name", "-pn", required=False, default=None, help="Project (bucket) name"
)
opt_project_req = click.option("--project-name", "-pn", required=True, help="Project (bucket) name")
opt_input_file = click.option("--input", "-i", "input_file", required=True, help="Path to input CSV")
opt_ins_input_file = click.option("--input", "-i", "input_file", required=False, help="Path to input CSV")
opt_output_file = click.option(
  "--output", "-o", "output_file", required=False, default=None, help="Path to output file (csv/h5)"
)
opt_access = click.option(
  "--access",
  type=click.Choice(["both", "public", "private"]),
  default=None,
  required=True,
  show_default=True,
  help="Which buckets to search when project-name not provided",
)
opt_yes_flag = click.option("--yes", "-y", is_flag=True, help="Confirm deletion")
opt_force_flag = click.option(
  "--force",
  is_flag=True,
  help="Allow writes to default isaura-public or isaura-private buckets",
)
opt_dump_outdir = click.option("--output-dir", "-o", required=False, help="Local output directory")
opt_approx = click.option(
  "--approximate",
  "-nn",
  is_flag=True,
  default=False,
  help="Use Approximate Nearest Neighbor search for result retrieval. [red bold]Under development — may return incomplete or unexpected results.[/]",
)
opt_cloud = click.option("--remote", "-r", "cloud", is_flag=True, default=False, help="Use the remote store instead of local")
opt_start = click.option("--start", "-s", is_flag=True, default=False, help="Start local services (MinIO)")
opt_stop = click.option("--stop", is_flag=True, default=False, help="Stop local services")
opt_status = click.option("--status", is_flag=True, default=False, help="Show status of local services")
opt_isaura_dir = click.option(
  "--isaura-dir",
  "-d",
  required=False,
  default=".",
  show_default=True,
  help="Path to an isaura folder (used mainly to resolve output defaults).",
)
opt_stats_outdir = click.option(
  "--output-dir", "-o", required=True, help="Folder where the stats JSON will be written."
)


@cli.command("configure")
@click.option("--remote", is_flag=True, default=False, help="Add or update remote/cloud credentials")
@click.option("--update", "do_update", is_flag=True, default=False, help="Update an existing credential interactively")
@click.option("--show-secrets", is_flag=True, default=False, help="Show all credential values unmasked")
@click.option("--test-credentials", "test_creds", is_flag=True, default=False, help="Test local and cloud credential connectivity")
def configure(remote, do_update, show_secrets, test_creds):
  """Show or update isaura configuration."""
  from isaura.configure import configure_show, configure_remote_interactive, configure_update_interactive, configure_test_credentials
  if test_creds:
    configure_test_credentials()
  elif do_update:
    configure_update_interactive()
  elif remote:
    configure_remote_interactive()
  else:
    configure_show(reveal=show_secrets)


@cli.command("create")
@click.option("--project-name", "-pn", required=True, help="Name for the new project bucket")
@click.option("--access", type=click.Choice(["public", "private"]), required=True, help="Access level for molecules stored in this bucket")
def create(project_name, access):
  """Create a new local project bucket."""
  import re
  import json
  import tempfile
  from isaura.base import MinioStore
  from isaura.const import (
    MINIO_ENDPOINT, MINIO_LOCAL_AK, MINIO_LOCAL_SK,
    DEFAULT_BUCKET_NAME, DEFAULT_PRIVATE_BUCKET_NAME, ACCESS_FILE,
  )
  if not re.fullmatch(r"[a-zA-Z0-9-]+", project_name):
    console.print("[red]Invalid project name.[/] Only alphanumeric characters and [bold]-[/bold] are allowed.")
    sys.exit(1)
  if project_name in (DEFAULT_BUCKET_NAME, DEFAULT_PRIVATE_BUCKET_NAME):
    console.print(f"[red]The name [bold]{project_name}[/bold] is reserved and cannot be used.[/]")
    sys.exit(1)
  store = MinioStore(endpoint=MINIO_ENDPOINT, access=MINIO_LOCAL_AK, secret=MINIO_LOCAL_SK)
  store.ensure_bucket(project_name)
  with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
    json.dump({"access": access}, f)
    tmp_path = f.name
  store.upload_file(tmp_path, project_name, ACCESS_FILE)
  console.print(f"[green]Project [bold]{project_name}[/bold] created with access=[bold]{access}[/bold].[/]")


@cli.command("write")
@apply_opts(opt_input_file, opt_project_req, opt_model, opt_version)
@click.option("--verbose", "-V", is_flag=True, default=False, help="Show detailed internal logs")
def write(input_file, project_name, model, version, verbose):
  """Store model outputs in a project bucket."""
  if verbose:
    logger.set_verbosity(True)
  _check_model_filename(model, input_file)
  with IsauraWriter(
    input_csv=input_file, model_id=model, model_version=version, bucket=project_name
  ) as w:
    w.write()


@cli.command("read")
@apply_opts(opt_input_file, opt_project_req, opt_model, opt_output_file)
@click.option("--version", "-v", required=False, default=None, help="Model version (default: latest stored version)")
# @apply_opts(opt_approx)  # approximate search disabled — under development
@click.option("--verbose", "-V", is_flag=True, default=False, help="Show detailed internal logs")
def read(input_file, project_name, model, output_file, version, verbose):
  """Retrieve stored model outputs for a set of inputs."""
  if verbose:
    logger.set_verbosity(True)
  _check_model_filename(model, output_file)
  if version is None:
    from isaura.base import MinioStore
    from isaura.const import MINIO_ENDPOINT, MINIO_LOCAL_AK, MINIO_LOCAL_SK
    try:
      _store = MinioStore(endpoint=MINIO_ENDPOINT, access=MINIO_LOCAL_AK, secret=MINIO_LOCAL_SK)
      version = _resolve_version(model, project_name, _store)
      console.print(f"[dim]No version specified — using latest: {version}[/dim]")
    except Exception:
      version = "v1"
  # approximate = False  # disabled — uncomment opt_approx above to re-enable
  with IsauraReader(
    model_id=model, model_version=version, bucket=project_name, input_csv=input_file, approximate=False
  ) as r:
    total = sum(len(chunk) for chunk in r.read_batched(output_csv=output_file))
  dest = f" → {output_file}" if output_file else ""
  console.print(f"[green]✓[/green] [bold]{model}/{version}[/bold]: [bold]{total}[/bold] rows{dest}")


@cli.command("pull")
@apply_opts(opt_input_file, opt_project_req, opt_model)
@click.option("--version", "-v", required=False, default=None, help="Model version (default: latest stored version in the remote bucket)")
@click.option("--verbose", "-V", is_flag=True, default=False, help="Show detailed internal logs")
def pull(input_file, project_name, model, version, verbose):
  """Pull model outputs from the cloud store to local."""
  if verbose:
    logger.set_verbosity(True)
  pn = project_name
  if version is None:
    from isaura.base import MinioStore
    from isaura.const import MINIO_ENDPOINT_CLOUD, MINIO_CLOUD_AK, MINIO_CLOUD_SK
    try:
      _store = MinioStore(endpoint=MINIO_ENDPOINT_CLOUD, access=MINIO_CLOUD_AK, secret=MINIO_CLOUD_SK)
      version = _resolve_version(model, pn, _store)
      console.print(f"[dim]No version specified — using latest: {version}[/dim]")
    except Exception:
      version = "v1"
  with IsauraPull(model_id=model, model_version=version, bucket=pn, input_csv=input_file) as pl:
    pl.pull()


@cli.command("push")
@apply_opts(opt_project, opt_model, opt_version)
@click.option("--verbose", "-V", is_flag=True, default=False, help="Show detailed internal logs")
def push(project_name, model, version, verbose):
  """Push local model outputs to the cloud store."""
  if verbose:
    logger.set_verbosity(True)
  p = IsauraPush(model, version, project_name)
  p.push()


@cli.command("persist")
@apply_opts(opt_model, opt_version, opt_project_req, opt_dump_outdir)
@click.option("--verbose", "-V", is_flag=True, default=False, help="Show detailed internal logs")
def cp(model, version, project_name, output_dir, verbose):
  """Copy model outputs from a project bucket into the canonical isaura-public/private buckets."""
  if verbose:
    logger.set_verbosity(True)
  with IsauraCopy(model_id=model, model_version=version, bucket=project_name, output_dir=output_dir) as c:
    priv, pub = c.copy()
  if output_dir is None:
    parts = []
    if pub:
      parts.append(f"[bold]{pub}[/bold] → isaura-public")
    if priv:
      parts.append(f"[bold]{priv}[/bold] → isaura-private")
    console.print(f"[green]✓[/green] [bold]{model}/{version}[/bold] copied: {', '.join(parts) if parts else 'nothing new'}")


# @cli.command("move")
# @apply_opts(opt_model, opt_version, opt_project_req)
# def mv(model, version, project_name):
#   """Move model outputs into the canonical isaura-public/private buckets (removes source)."""
#   with IsauraMover(model_id=model, model_version=version, bucket=project_name) as m:
#     m.move()
#     logger.info(f"Move done for {model}/{version} from {project_name}")


@cli.command("engine")
@apply_opts(opt_start, opt_stop, opt_status)
def engine(start, stop, status):
  """Manage local isaura services (MinIO). Use --start, --stop, or --status."""
  if stop:
    run_docker_compose(up=False)
  elif start:
    run_docker_compose(up=True)
  else:
    info = get_engine_status()
    table = make_table("Local services", [
      {"key": "name", "name": "Service", "justify": "left", "style": "bold"},
      {"key": "status", "name": "Status", "justify": "left"},
    ], [{"name": "Docker", "status": "[green]running[/]" if info["docker"] else "[red]not running[/]"}] + [
      {"name": c["name"], "status": c["status"]} for c in info["containers"]
    ])
    console.print(table)
    if not info["docker"]:
      console.print("[dim]Open Docker Desktop and run [bold]isaura engine --start[/bold] to start local services.[/]")


@cli.command("remove")
@click.option("--model-id", "-m", "model", required=True, help="Ersilia model ID (eosxxxx)")
@click.option("--version", "-v", default="v1", show_default=True, help="Model version")
@click.option("--project-name", "-pn", required=True, help="Project (bucket) name")
@click.option("--input", "-i", "input_file", required=True, help="CSV of molecules to remove")
@click.option("--verbose", "-V", is_flag=True, default=False, help="Show detailed internal logs")
@click.option("--yes", "-y", is_flag=True, default=False, help="Skip confirmation prompt")
def rm(model, version, project_name, input_file, verbose, yes):
  """Remove specific molecules from a project bucket."""
  if verbose:
    logger.set_verbosity(True)
  if not yes:
    console.print(
      f"[yellow]Warning:[/yellow] This will permanently remove molecules from "
      f"[bold]{model}/{version}[/bold] in [bold]{project_name}[/bold] and cannot be undone."
    )
    if not click.confirm("Proceed?", default=False):
      console.print("[dim]Cancelled.[/dim]")
      return
  with IsauraMolRemover(
    model_id=model, model_version=version, bucket=project_name, input_csv=input_file
  ) as r:
    n_removed, n_not_found = r.remove()
  not_found_str = f", [yellow]{n_not_found} not found[/yellow]" if n_not_found else ""
  console.print(
    f"[green]✓[/green] [bold]{model}/{version}[/bold] → [bold]{project_name}[/bold]: "
    f"[bold]{n_removed}[/bold] molecules removed{not_found_str}"
  )


@cli.command("destroy")
@click.option("--project-name", "-pn", required=True, help="Project (bucket) name to destroy")
@click.option("--model-id", "-m", "model", required=False, default=None, help="If set, destroy only this model's data. Requires --version.")
@click.option("--version", "-v", required=False, default=None, show_default=True, help="Model version to destroy. Required when --model-id is set.")
@click.option("--yes", "-y", is_flag=True, default=False, help="Skip confirmation prompt")
def destroy(project_name, model, version, yes):
  """Destroy a project bucket or a specific model version within it.

  \b
  • Without --model-id: destroys the entire project bucket.
  • With --model-id and --version: destroys only that model/version combination.

  Cannot target isaura-public or isaura-private.
  """
  if model and not version:
    console.print("[red]--version is required when --model-id is set.[/] Example: [bold]-m eos4e40 -v v1[/bold]")
    sys.exit(1)
  from isaura.base import MinioStore
  from isaura.const import MINIO_ENDPOINT, MINIO_LOCAL_AK, MINIO_LOCAL_SK
  if not model and project_name in (DEFAULT_BUCKET_NAME, DEFAULT_PRIVATE_BUCKET_NAME):
    console.print(f"[red][bold]{project_name}[/bold] is a reserved bucket and cannot be destroyed.[/] Use [bold]--model-id[/bold] and [bold]--version[/bold] to remove a specific model.")
    sys.exit(1)
  store = MinioStore(endpoint=MINIO_ENDPOINT, access=MINIO_LOCAL_AK, secret=MINIO_LOCAL_SK)
  store.require_bucket(project_name)
  if not yes:
    if model:
      console.print(
        f"[yellow]Warning:[/yellow] This will permanently delete all data for "
        f"[bold]{model}/{version}[/bold] in [bold]{project_name}[/bold] and cannot be undone."
      )
    else:
      console.print(
        f"[yellow]Warning:[/yellow] This will permanently destroy the entire project "
        f"[bold]{project_name}[/bold] and all its data and cannot be undone."
      )
    if not click.confirm("Proceed?", default=False):
      console.print("[dim]Cancelled.[/dim]")
      return
  if model:
    with console.status(f"Deleting {model}/{version} from {project_name}...", spinner="dots"):
      n = store.delete_prefix(project_name, f"{model}/{version}/")
    console.print(f"[green]✓[/green] [bold]{model}/{version}[/bold] deleted from [bold]{project_name}[/bold]: {n} objects removed")
  else:
    with console.status(f"Destroying {project_name}...", spinner="dots"):
      n = store.delete_prefix(project_name, "")
      store.client.delete_bucket(Bucket=project_name)
    console.print(f"[green]✓[/green] Project [bold]{project_name}[/bold] destroyed: {n} objects removed")


@cli.command("inspect")
@click.option("--model-id", "-m", "model_id", required=True, help="Ersilia model ID to inspect (e.g. eos4e40)")
@click.option("--version", "-v", default=None, help="Model version (default: latest stored version)")
@click.option("--project-name", "-pn", required=True, help="Project bucket to search (e.g. isaura-public, isaura-private).")
@click.option("--input", "-i", "input_file", required=True,
              help="CSV of molecules to check against the store.")
@click.option("--output", "-o", "output_file", required=True,
              help="Path to write the results CSV.")
@click.option("--remote", "-r", "cloud", is_flag=True, default=False,
              help="Query the remote (cloud) store. Default: local MinIO instance.")
def cmd_inspect(model_id, version, project_name, input_file, output_file, cloud):
  """Check which molecules from an input CSV are already stored for a model."""
  from isaura.base import MinioStore
  if cloud:
    from isaura.const import MINIO_ENDPOINT_CLOUD, MINIO_CLOUD_AK, MINIO_CLOUD_SK
    _store = MinioStore(endpoint=MINIO_ENDPOINT_CLOUD, access=MINIO_CLOUD_AK, secret=MINIO_CLOUD_SK)
  else:
    from isaura.const import MINIO_ENDPOINT, MINIO_LOCAL_AK, MINIO_LOCAL_SK
    _store = MinioStore(endpoint=MINIO_ENDPOINT, access=MINIO_LOCAL_AK, secret=MINIO_LOCAL_SK)
  _store.require_bucket(project_name)
  if version is None:
    try:
      version = _resolve_version(model_id, project_name, _store)
      console.print(f"[dim]No version specified — using latest: {version}[/dim]")
    except Exception:
      version = "v1"
  insp = IsauraInspect(
    model_id=model_id, model_version=version, project_name=project_name, cloud=cloud
  )
  if cloud:
    for b in insp.buckets():
      insp._clients(b)
    ctx = console.status("Inspecting...", spinner="dots")
  else:
    ctx = contextlib.nullcontext()
  with ctx:
    df = insp.inspect_inputs(input_file, output_file)
  logger.info(f"wrote {len(df)} rows → {output_file}")


@cli.command("catalog")
@click.option("--remote", "-r", "cloud", is_flag=True, default=False,
              help="List models in the remote store instead of local")
@click.option("--project-name", "-pn", required=False, default=None,
              help="Project (bucket) name — e.g. isaura-public or isaura-private")
def cmd_inspect_models(project_name, cloud):
  """List all models stored in a project bucket."""
  if not project_name:
    if cloud:
      console.print(
        "[yellow]Please specify a project name.[/] "
        "Use [bold]-pn isaura-public[/bold] or [bold]-pn isaura-private[/bold]."
      )
    else:
      console.print(
        "[yellow]Please specify a project name.[/] "
        "Run [bold]isaura info[/bold] to see available local projects, "
        "then re-run with [bold]-pn <project>[/bold]."
      )
    return
  if cloud:
    from isaura.const import MINIO_CLOUD_AK
    if not MINIO_CLOUD_AK:
      console.print(
        "[yellow]No remote credentials configured.[/] "
        "Run [bold]isaura configure --remote[/bold] first."
      )
      return
  insp = IsauraInspect(cloud=cloud)
  insp._clients(project_name)
  try:
    with console.status("Fetching catalog...", spinner="dots"):
      rows = insp.inspect_models(project_name, prefix_filter="")
  except Exception as e:
    msg = str(e)
    if "NoSuchBucket" in msg or "does not exist" in msg.lower():
      console.print(f"[red]Project [bold]{project_name}[/bold] not found.[/] Run [bold]isaura info{'  --remote' if cloud else ''}[/bold] to see available projects.")
    elif "InvalidAccessKeyId" in msg or "SignatureDoesNotMatch" in msg or "Access Key" in msg:
      console.print("[red]Invalid credentials.[/] Run [bold]isaura configure --remote[/bold] to update them.")
    elif "connect" in msg.lower() or "endpoint" in msg.lower():
      target = "remote store" if cloud else "local MinIO"
      console.print(f"[red]Could not connect to {target}.[/]" + ("" if cloud else " Is Docker running? Try [bold]isaura engine --start[/bold]."))
    else:
      console.print(f"[red]Catalog fetch failed:[/] {msg}")
    return
  if not rows:
    console.print(f"[yellow]No models found in {project_name}[/]")
    return
  cols = inspect_table_cloud if cloud else inspect_table
  table = make_table(f"Model catalog in {project_name}", cols, rows)
  console.print(table)


@cli.command("info")
@click.option("--remote", "-r", "cloud", is_flag=True, default=False, help="List projects in the remote store instead of local")
def info(cloud):
  """List available projects (buckets) in the local or remote store."""
  from isaura.base import MinioStore
  from isaura.const import (
    MINIO_ENDPOINT, MINIO_LOCAL_AK, MINIO_LOCAL_SK,
    MINIO_ENDPOINT_CLOUD, MINIO_CLOUD_AK, MINIO_CLOUD_SK,
  )
  if cloud:
    if not MINIO_CLOUD_AK:
      console.print("[yellow]No cloud credentials configured.[/] Run [bold]isaura configure --remote[/bold] first.")
      return
    endpoint, access, secret = MINIO_ENDPOINT_CLOUD, MINIO_CLOUD_AK, MINIO_CLOUD_SK
    title = "Remote projects"
  else:
    endpoint, access, secret = MINIO_ENDPOINT, MINIO_LOCAL_AK, MINIO_LOCAL_SK
    title = "Local projects"
  try:
    store = MinioStore(endpoint=endpoint, access=access, secret=secret)
  except SystemExit:
    target = "remote store" if cloud else "local MinIO"
    console.print(f"[red]Could not connect to {target}.[/]" + ("" if cloud else " Is Docker running? Try [bold]isaura engine --start[/bold]."))
    return
  buckets = store.client.list_buckets().get("Buckets", [])
  if not buckets:
    console.print(f"[yellow]No projects found in {title.lower()}.[/]")
    return

  def _access_label(bucket_name):
    import json, tempfile
    if "public" in bucket_name:
      return "[green]public[/green]"
    if "private" in bucket_name:
      return "[red]private[/red]"
    try:
      tmp = tempfile.mktemp()
      store.download_file(bucket_name, "access.json", tmp)
      with open(tmp) as f:
        val = json.load(f).get("access", "")
      if val == "public":
        return "[green]public[/green]"
      if val == "private":
        return "[red]private[/red]"
    except Exception:
      pass
    return "[dim]–[/dim]"

  def _location(bucket_name):
    if cloud:
      return f"{endpoint}/{bucket_name}"
    return os.path.expanduser(f"~/minio-data/{bucket_name}")

  rows = [{"name": b["Name"], "created": b["CreationDate"].strftime("%Y-%m-%d %H:%M"), "access": _access_label(b["Name"]), "location": _location(b["Name"])} for b in buckets]
  table = make_table(title, [
    {"key": "name", "name": "Project", "justify": "left", "style": "bold"},
    {"key": "access", "name": "Access", "justify": "left"},
    {"key": "created", "name": "Created", "justify": "left"},
    {"key": "location", "name": "URL" if cloud else "Path", "justify": "left", "style": "dim"},
  ], rows)
  console.print(table)


@cli.command("stats")
@click.option("--project-name", "-pn", required=True,
              help="Project (bucket) name to collect statistics for (e.g. isaura-public, isaura-private, or a custom project).")
@click.option("--remote", "-r", "cloud", is_flag=True, default=False,
              help="Collect statistics from the remote (cloud) store. Default: local MinIO instance.")
@click.option("--output-dir", "-o", required=True,
              help="Folder where the stats JSON file will be written. "
                   "The file is named isaura_stats_<timestamp>.json.")
@click.option("--isaura-dir", "-d", required=False, default=".", hidden=True,
              help="Path to an isaura folder (used mainly to resolve output defaults).")
def cmd_stats(project_name, cloud, output_dir, isaura_dir):
  """Generate a JSON inventory of all models stored in a project bucket."""
  from isaura.base import MinioStore
  if cloud:
    from isaura.const import MINIO_ENDPOINT_CLOUD, MINIO_CLOUD_AK, MINIO_CLOUD_SK
    MinioStore(endpoint=MINIO_ENDPOINT_CLOUD, access=MINIO_CLOUD_AK, secret=MINIO_CLOUD_SK).require_bucket(project_name)
  else:
    from isaura.const import MINIO_ENDPOINT, MINIO_LOCAL_AK, MINIO_LOCAL_SK
    MinioStore(endpoint=MINIO_ENDPOINT, access=MINIO_LOCAL_AK, secret=MINIO_LOCAL_SK).require_bucket(project_name)
  ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
  out_path = os.path.join(output_dir, f"isaura_stats_{ts}.json")
  st = IsauraStat(project_name=project_name, cloud=cloud, endpoint=None)
  written = st.write_json(out_path)
  logger.info(f"isaura stats wrote: {written}")


if __name__ == "__main__":
  cli()
