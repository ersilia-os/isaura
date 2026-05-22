import contextlib, datetime, os, sys
import rich_click as click
import rich_click.rich_click as rc
from isaura.manage import (
  IsauraMover,
  IsauraCopy,
  IsauraRemover,
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
  spinner,
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


@click.group()
def cli():
  pass


show_figlet()
opt_model = click.option("--model", "-m", required=True, help="Ersilia model id (eosxxxx)")
opt_model_opt = click.option("--model", "-m", required=False, default=None, help="Ersilia model id (eosxxxx)")
opt_version = click.option("--version", "-v", default="v1", show_default=True, help="Model version")
opt_project = click.option(
  "--project-name", "-pn", required=False, default=None, help="Project (bucket) name"
)
opt_project_req = click.option("--project-name", "-pn", required=True, help="Project (bucket) name")
opt_input_file = click.option("--input-file", "-i", required=True, help="Path to input CSV")
opt_ins_input_file = click.option("--input-file", "-i", required=False, help="Path to input CSV")
opt_output_file = click.option(
  "--output-file", "-o", required=False, default=None, help="Path to output file (csv/h5)"
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


@cli.command("write")
@apply_opts(opt_input_file, opt_project, opt_access, opt_model, opt_version, opt_force_flag)
def write(input_file, project_name, access, model, version, force):
  """Store model outputs in a project bucket."""
  if project_name in (DEFAULT_PRIVATE_BUCKET_NAME, DEFAULT_BUCKET_NAME):
    if not force:
      logger.error("Access denied to write to default project names. Re-run with --force to allow it.")
      sys.exit(1)
    logger.warning(f"Force-enabled write directly to protected bucket: {project_name}")
  if project_name is None:
    logger.error("Please specify the project name in order to write the data!")
    sys.exit(1)
  with IsauraWriter(
    input_csv=input_file, model_id=model, model_version=version, bucket=project_name, access=access
  ) as w:
    w.write()


@cli.command("read")
@apply_opts(opt_input_file, opt_project, opt_model, opt_version, opt_output_file, opt_approx)
def read(input_file, project_name, model, version, output_file, approximate):
  """Retrieve stored model outputs for a set of inputs."""
  if approximate:
    logger.warning("Approximate Nearest Neighbor search is under active development and may return incomplete or unexpected results.")
  with IsauraReader(
    model_id=model, model_version=version, bucket=project_name, input_csv=input_file, approximate=approximate
  ) as r:
    for _ in r.read_batched(output_csv=output_file):
      pass


@cli.command("pull")
@apply_opts(opt_input_file, opt_project, opt_model, opt_version)
def pull(input_file, project_name, model, version):
  """Pull model outputs from the cloud store to local."""
  pn = project_name or DEFAULT_BUCKET_NAME
  with IsauraPull(model_id=model, model_version=version, bucket=pn, input_csv=input_file) as pl:
    pl.pull()


@cli.command("push")
@apply_opts(opt_project, opt_model, opt_version)
def push(project_name, model, version):
  """Push local model outputs to the cloud store."""
  p = IsauraPush(model, version, project_name)
  p.push()


@cli.command("copy")
@apply_opts(opt_model, opt_version, opt_project_req, opt_dump_outdir)
def cp(model, version, project_name, output_dir):
  """Copy model outputs from a project bucket into the canonical isaura-public/private buckets."""
  with IsauraCopy(model_id=model, model_version=version, bucket=project_name, output_dir=output_dir) as c:
    if output_dir is None:
      priv, pub = c.copy()
      logger.info(f"Copied private_new={priv} public_new={pub} from {project_name}")
    else:
      c.copy()


@cli.command("move")
@apply_opts(opt_model, opt_version, opt_project_req)
def mv(model, version, project_name):
  """Move model outputs into the canonical isaura-public/private buckets (removes source)."""
  with IsauraMover(model_id=model, model_version=version, bucket=project_name) as m:
    m.move()
    logger.info(f"Move done for {model}/{version} from {project_name}")


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
@apply_opts(opt_model_opt, opt_version, opt_project_req, opt_yes_flag)
def rm(model, version, project_name, yes):
  """Delete model outputs from a bucket. Requires --yes to confirm."""
  if not yes:
    logger.info("Add --yes to confirm deletion")
    sys.exit(1)
  if model:
    with IsauraRemover(model_id=model, model_version=version, bucket=project_name) as r:
      r.remove()
    logger.info(f"Remove done for {model}/{version} in {project_name}")
  else:
    with IsauraRemover(project_name=project_name) as r:
      r.remove()
    logger.info(f"Remove done for all data in {project_name}")


@cli.command("inspect")
@apply_opts(opt_model, opt_version, opt_project, opt_access, opt_ins_input_file, opt_output_file, opt_cloud)
@click.argument("what", type=click.Choice(["inputs"]), required=False, default="inputs")
def cmd_inspect(what, model, version, project_name, access, input_file, output_file, cloud):
  """List available inputs or entries for a model in the store."""
  insp = IsauraInspect(
    model_id=model, model_version=version, project_name=project_name, access=access, cloud=cloud
  )
  if cloud:
    for b in insp.buckets():
      insp._clients(b)
    ctx = console.status("Inspecting...", spinner="dots")
  else:
    ctx = contextlib.nullcontext()
  with ctx:
    if input_file:
      df = insp.inspect_inputs(input_file, output_file)
    else:
      df = insp.list_available(output_file)
  logger.info(f"wrote {len(df)} rows{(' -> ' + output_file if output_file else '')}")


@cli.command("catalog")
@apply_opts(opt_project, opt_cloud)
def cmd_inspect_models(project_name, cloud):
  """List all models stored in a project bucket."""
  insp = IsauraInspect(cloud=cloud)
  # warm up the client so the health-check log fires before the spinner
  insp._clients(project_name)
  with console.status("Fetching catalog...", spinner="dots"):
    rows = insp.inspect_models(project_name, prefix_filter="")
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
  rows = [{"name": b["Name"], "created": b["CreationDate"].strftime("%Y-%m-%d %H:%M")} for b in buckets]
  table = make_table(title, [
    {"key": "name", "name": "Project", "justify": "left", "style": "bold"},
    {"key": "created", "name": "Created", "justify": "left"},
  ], rows)
  console.print(table)


@cli.command("stats")
@apply_opts(opt_project, opt_access, opt_cloud, opt_stats_outdir, opt_isaura_dir)
def cmd_stats(project_name, access, cloud, output_dir, isaura_dir):
  """Generate a JSON inventory of all stored models and their sizes."""
  ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
  out_path = os.path.join(output_dir, f"isaura_stats_{ts}.json")
  st = IsauraStat(project_name=project_name, access=access, cloud=cloud, endpoint=None)
  written = st.write_json(out_path)
  logger.info(f"isaura stats wrote: {written}")


if __name__ == "__main__":
  cli()
