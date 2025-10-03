import sys

import rich_click as click
import rich_click.rich_click as rc

from isaura.manage import (
  IsauraMover,
  IsauraCopy,
  IsauraRemover,
  IsauraCleaner,
  IsauraWriter,
  IsauraReader,
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


opt_model = click.option("--model", "-m", required=True, help="Ersilia model id (eosxxxx)")
opt_version = click.option("--version", "-v", default="v1", show_default=True, help="Model version")
opt_project_opt = click.option(
  "--project-name", "-pn", required=False, default=None, help="Project (bucket) name"
)
opt_project_req = click.option("--project-name", "-pn", required=True, help="Project (bucket) name")
opt_access = click.option(
  "--access-level", "-al", default="public", show_default=True, help="public | private"
)
opt_input_file = click.option("--input-file", "-i", required=True, help="Path to input CSV")
opt_output_file = click.option(
  "--output-file",
  "-o",
  required=False,
  default=None,
  help="Path to output file (csv/h5)",
)
opt_metadata = click.option(
  "--metadata", "-md", required=False, default=None, help="Path to metadata.json"
)
opt_yes_flag = click.option("--yes", "-y", is_flag=True, help="Confirm deletion")
opt_dump_outdir = click.option("--output-dir", "-o", required=True, help="Local output directory")


def apply_opts(*opts):
  def _wrap(f):
    for opt in reversed(opts):
      f = opt(f)
    return f

  return _wrap


@click.group()
def cli():
  pass


@cli.command("write")
@apply_opts(opt_input_file, opt_project_opt, opt_access, opt_model, opt_version, opt_metadata)
def write(input_file, project_name, access_level, model, version, metadata):
  with IsauraWriter(
    input_csv=input_file,
    model_id=model,
    model_version=version,
    bucket=project_name,
    acess_level=access_level,
    metadata_path=metadata,
  ) as w:
    w.write()


@cli.command("read")
@apply_opts(opt_input_file, opt_project_opt, opt_access, opt_model, opt_version, opt_output_file)
def read(input_file, project_name, access_level, model, version, output_file):
  r = IsauraReader(model_id=model, model_version=version, bucket=project_name, input_csv=input_file)
  r.read(output_csv=output_file)


@cli.command("cp")
@apply_opts(opt_model, opt_version, opt_project_req)
def cp(model, version, project_name):
  c = IsauraCopy(model_id=model, model_version=version, project_name=project_name)
  priv, pub = c.copy()
  click.echo(f"Copied private_new={priv} public_new={pub} from {project_name}")


@cli.command("mv")
@apply_opts(opt_model, opt_version, opt_project_req)
def mv(model, version, project_name):
  m = IsauraMover(model_id=model, model_version=version, project_name=project_name)
  m.move()
  click.echo(f"Move done for {model}/{version} from {project_name}")


@cli.command("rm")
@apply_opts(opt_model, opt_version, opt_project_req, opt_yes_flag)
def rm(model, version, project_name, yes):
  if not yes:
    click.echo("Add --yes to confirm deletion")
    sys.exit(1)
  r = IsauraRemover(model_id=model, model_version=version, project_name=project_name)
  r.remove()
  click.echo(f"Remove done for {model}/{version} in {project_name}")


@cli.command("clean")
@apply_opts(opt_model, opt_version)
def clean(model, version):
  c = IsauraCleaner(model_id=model, model_version=version)
  res = c.clean()
  for bucket, n in res.items():
    click.echo(f"Cleaned {bucket}:{model}/{version}/tranches -> deleted {n} objects")


if __name__ == "__main__":
  cli()
