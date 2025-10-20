import sys

import rich_click as click
import rich_click.rich_click as rc

from isaura.manage import (
  IsauraMover,
  IsauraCopy,
  IsauraRemover,
  IsauraWriter,
  IsauraReader,
  IsauraInspect,
)
from isaura.helpers import logger, console, make_table, inspect_table

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


opt_model = click.option("--model", "-m", required=True, help="Ersilia model id (eosxxxx)")
opt_version = click.option("--version", "-v", default="v1", show_default=True, help="Model version")
opt_project = click.option(
  "--project-name", "-pn", required=False, default=None, help="Project (bucket) name"
)
opt_project_req = click.option("--project-name", "-pn", required=True, help="Project (bucket) name")
opt_input_file = click.option("--input-file", "-i", required=True, help="Path to input CSV")
opt_output_file = click.option(
  "--output-file",
  "-o",
  required=False,
  default=None,
  help="Path to output file (csv/h5)",
)
opt_access = click.option(
  "--access",
  type=click.Choice(["public", "private", "both"]),
  default=None,
  show_default=True,
  help="Which buckets to search when project-name not provided",
)
opt_yes_flag = click.option("--yes", "-y", is_flag=True, help="Confirm deletion")
opt_dump_outdir = click.option("--output-dir", "-o", required=False, help="Local output directory")
opt_approx = click.option(
  "--approximate",
  "-a",
  is_flag=True,
  default=False,
  help="Specifies whether to use Approximate Nearest Neighbor search for result retrieval or not.",
)


@cli.command("write")
@apply_opts(opt_input_file, opt_project, opt_access, opt_model, opt_version)
def write(input_file, project_name, access, model, version):
  with IsauraWriter(
    input_csv=input_file,
    model_id=model,
    model_version=version,
    bucket=project_name,
    access=access,
  ) as w:
    w.write()


@cli.command("read")
@apply_opts(opt_input_file, opt_project, opt_model, opt_version, opt_output_file, opt_approx)
def read(input_file, project_name, model, version, output_file, approximate):
  r = IsauraReader(
    model_id=model, model_version=version, bucket=project_name, input_csv=input_file, approximate=approximate
  )
  r.read(output_csv=output_file)


@cli.command("copy")
@apply_opts(opt_model, opt_version, opt_project_req, opt_dump_outdir)
def cp(model, version, project_name, output_dir):
  c = IsauraCopy(model_id=model, model_version=version, bucket=project_name, output_dir=output_dir)
  if output_dir is None:
    priv, pub = c.copy()
    logger.info(f"Copied private_new={priv} public_new={pub} from {project_name}")
  else:
    c.copy()


@cli.command("move")
@apply_opts(opt_model, opt_version, opt_project_req)
def mv(model, version, project_name):
  m = IsauraMover(model_id=model, model_version=version, project_name=project_name)
  m.move()
  logger.info(f"Move done for {model}/{version} from {project_name}")


@cli.command("remove")
@apply_opts(opt_model, opt_version, opt_project_req, opt_yes_flag)
def rm(model, version, project_name, yes):
  if not yes:
    logger.info("Add --yes to confirm deletion")
    sys.exit(1)
  r = IsauraRemover(model_id=model, model_version=version, project_name=project_name)
  r.remove()
  logger.info(f"Remove done for {model}/{version} in {project_name}")


@cli.command("inspect")
@apply_opts(opt_model, opt_version, opt_project, opt_access, opt_input_file, opt_output_file)
@click.argument("what", type=click.Choice(["inputs"]), required=False, default="inputs")
def cmd_inspect(what, model, version, project_name, access, input_file, output_file):
  insp = IsauraInspect(model_id=model, model_version=version, project_name=project_name, access=access)
  if input_file:
    df = insp.inspect_inputs(input_file, output_file)
  else:
    df = insp.list_available(output_file)
  logger.info(f"wrote {len(df)} rows{(' -> ' + output_file) if output_file else ''}")


@cli.command("catalog")
@apply_opts(opt_project)
@click.option("-f", "--filter", default="", required=False, help="Model id prefix filter")
def cmd_inspect_models(project_name, filter):
  insp = IsauraInspect(model_id="_", model_version="_")
  rows = insp.inspect_models(project_name, prefix_filter=filter)
  if not rows:
    console.print(f"[yellow]No models found in {project_name}[/]")
    return
  table = make_table(
    f"Models in {project_name}",
    inspect_table,
    rows,
  )
  console.print(table)


if __name__ == "__main__":
  cli()
