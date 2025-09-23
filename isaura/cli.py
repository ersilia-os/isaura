import rich_click as click
import rich_click.rich_click as rc
from isaura.register.local import IsauraWriter

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


@click.group()
def cli():
  pass


@cli.command()
@click.option("--input-file", "-i", required=True, help="Path to the input file.")
@click.option(
  "--project-name",
  "-pn",
  required=False,
  default=None,
  help="A custom project name to create and store the precalculation.",
)
@click.option(
  "--access-level",
  "-al",
  required=False,
  default="public",
  help="Two types of access levels are allowed. One is the public to store public data and the other is private to store private data. The private access level requires a secrete key in the environment variable as [SECERETE]",
)
@click.option("--model", "-m", required=True, help="Ersilia model id (eosxxxx)")
@click.option("--version", "-v", required=False, default="v1", help="Mode version")
def write(input_file, project_name, access_level, model, version):
  with IsauraWriter(
    input_csv=input_file,
    model_id=model,
    model_version=version,
    bucket=project_name,
    acess_level=access_level,
  ) as w:
    w.write()


if __name__ == "__main__":
  cli()
