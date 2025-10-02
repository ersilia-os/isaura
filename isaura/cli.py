import rich_click as click
import rich_click.rich_click as rc
from isaura.register.local import IsauraWriter, IsauraReader

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
@click.option(
  "--input-file",
  "-i",
  required=True,
  help="Path to the input file. The input file should contains an Ersilia formatted output result [key, input, col1, col2,..].",
)
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
@click.option(
  "--metadata",
  "-md",
  required=False,
  default=None,
  help="This file contains the information about the access level of input calculation",
)
@click.option("--model", "-m", required=True, help="Ersilia model id (eosxxxx)")
@click.option("--version", "-v", required=False, default="v1", help="Mode version")
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


@cli.command()
@click.option(
  "--input-file",
  "-i",
  required=True,
  help="Path to the input file. This should contains a list of smiles.",
)
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
@click.option(
  "--output-file",
  "-o",
  required=False,
  default=None,
  help="The path to output file (csv or h5)",
)
def read(input_file, project_name, access_level, model, version, output_file):
  reader = IsauraReader(
    model_id=model, model_version=version, bucket=project_name, input_csv=input_file
  )
  reader.read(output_file=output_file)


if __name__ == "__main__":
  cli()
