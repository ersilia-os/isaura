"""Configuration helpers for the `isaura configure` CLI command."""

import sys
from pathlib import Path

import click
import questionary
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from isaura.const import (
    DEFAULT_BUCKET_NAME,
    DEFAULT_PRIVATE_BUCKET_NAME,
    MINIO_ENDPOINT,
    MINIO_LOCAL_AK,
    MINIO_LOCAL_SK,
    MINIO_ENDPOINT_CLOUD,
    MINIO_CLOUD_AK,
    MINIO_CLOUD_SK,
    MINIO_PRIV_CLOUD_AK,
    MINIO_PRIV_CLOUD_SK,
)

console = Console()

CONFIG_DIR = Path.home() / ".isaura"
ENV_PATH = CONFIG_DIR / ".env"

LOCAL_DEFAULTS = {
    "MINIO_ENDPOINT": "http://127.0.0.1:9000",
    "NNS_ENDPOINT": "http://127.0.0.1:8080",
    "DEFAULT_BUCKET_NAME": "isaura-public",
    "DEFAULT_PRIVATE_BUCKET_NAME": "isaura-private",
    "MINIO_LOCAL_AK": "minioadmin123",
    "MINIO_LOCAL_SK": "minioadmin1234",
}

_SECRET_TOKENS = {"AK", "SK", "KEY", "PASSWORD", "SECRET"}


def _is_secret(key: str) -> bool:
    upper = key.upper()
    return "CLOUD" in upper and any(tok in upper for tok in _SECRET_TOKENS)


def _mask(value: str) -> str:
    return "••••••••"


def _read_env_file(path: Path) -> dict:
    """Return key/value pairs from a .env file, ignoring comments and blanks."""
    result = {}
    if not path.exists():
        return result
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            k, _, v = line.partition("=")
            result[k.strip()] = v.strip().strip('"').strip("'")
    return result


def _write_env_file(path: Path, vars_dict: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{k}={v}" for k, v in vars_dict.items()]
    path.write_text("\n".join(lines) + "\n")


def configure_show(reveal: bool = False) -> None:
    """Display the current configuration as a table, masking secrets unless reveal=True."""
    if not ENV_PATH.exists():
        console.print(
            "[yellow]No configuration file found.[/] "
            "It will be created automatically the next time you run an isaura command."
        )
        return

    existing = _read_env_file(ENV_PATH)
    table = Table(title=f"isaura configuration  ({ENV_PATH})", show_lines=False)
    table.add_column("Variable", style="bold")
    table.add_column("Value")

    for key, value in existing.items():
        display = value if (reveal or not _is_secret(key)) else _mask(value)
        table.add_row(key, display)

    console.print(table)


def configure_remote_interactive() -> None:
    """Interactively prompt for remote/cloud credentials and merge into ~/.isaura/.env."""
    console.print(Panel(
        "You will be asked for your remote MinIO credentials.\n"
        "Secret keys are hidden as you type.\n"
        "Press [bold]Enter[/bold] to keep the current value where shown.",
        title="Remote configuration",
        border_style="blue",
    ))

    existing = _read_env_file(ENV_PATH)

    endpoint = click.prompt(
        "Cloud endpoint",
        default=existing.get("MINIO_ENDPOINT_CLOUD", "http://83.48.73.209:8080"),
    )
    pub_ak = click.prompt("Public access key", hide_input=True)
    pub_sk = click.prompt("Public secret key", hide_input=True)

    console.print("[dim]Private bucket credentials are optional — press Enter to skip.[/]")
    priv_ak = click.prompt("Private access key (optional)", default="", hide_input=True, show_default=False)
    priv_sk = click.prompt("Private secret key (optional)", default="", hide_input=True, show_default=False)

    updates = {
        "MINIO_ENDPOINT_CLOUD": endpoint,
        "MINIO_CLOUD_AK": pub_ak,
        "MINIO_CLOUD_SK": pub_sk,
    }
    if priv_ak:
        updates["MINIO_PRIV_CLOUD_AK"] = priv_ak
    if priv_sk:
        updates["MINIO_PRIV_CLOUD_SK"] = priv_sk

    merged = {**existing, **updates}
    _write_env_file(ENV_PATH, merged)
    console.print(f"[green]Remote credentials saved to {ENV_PATH}[/]")


def _try_credential_check(endpoint, access, secret) -> tuple[bool | None, str]:
    """Attempt a MinioStore credential check, returning (ok, message).

    Returns (None, message) if the server is unreachable (ping → sys.exit),
    (True, message) on success, (False, message) on auth failure.
    """
    from isaura.base import MinioStore
    try:
        store = MinioStore(endpoint=endpoint, access=access, secret=secret)
    except SystemExit:
        return None, "[red]✗ server unreachable[/]"
    except Exception:
        return None, "[red]✗ server unreachable[/]"
    ok = store.credential_check()
    return ok, ("[green]✓ connected[/]" if ok else "[red]✗ auth failed[/]")


def configure_test_credentials() -> None:
    """Test local and cloud MinIO credentials and print a result table."""
    from isaura.metadata import docker_is_running

    rows = []

    # Local
    if not docker_is_running():
        rows.append({"target": "Local", "bucket": DEFAULT_BUCKET_NAME, "result": "[red]✗ Docker not running[/]"})
    else:
        _, msg = _try_credential_check(MINIO_ENDPOINT, MINIO_LOCAL_AK, MINIO_LOCAL_SK)
        rows.append({"target": "Local", "bucket": DEFAULT_BUCKET_NAME, "result": msg})

    # Cloud public
    if not MINIO_CLOUD_AK:
        rows.append({"target": "Cloud public", "bucket": DEFAULT_BUCKET_NAME,
                     "result": "[yellow]✗ not configured[/]  run [bold]isaura configure --remote[/bold]"})
    else:
        _, msg = _try_credential_check(MINIO_ENDPOINT_CLOUD, MINIO_CLOUD_AK, MINIO_CLOUD_SK)
        rows.append({"target": "Cloud public", "bucket": DEFAULT_BUCKET_NAME, "result": msg})

    # Cloud private
    if not MINIO_PRIV_CLOUD_AK:
        rows.append({"target": "Cloud private", "bucket": DEFAULT_PRIVATE_BUCKET_NAME,
                     "result": "[yellow]✗ not configured[/]  run [bold]isaura configure --remote[/bold]"})
    else:
        _, msg = _try_credential_check(MINIO_ENDPOINT_CLOUD, MINIO_PRIV_CLOUD_AK, MINIO_PRIV_CLOUD_SK)
        rows.append({"target": "Cloud private", "bucket": DEFAULT_PRIVATE_BUCKET_NAME, "result": msg})

    table = Table(title="Credential check", show_lines=False)
    table.add_column("Target", style="bold")
    table.add_column("Bucket")
    table.add_column("Result")
    for row in rows:
        table.add_row(row["target"], row["bucket"], row["result"])
    console.print(table)


def configure_update_interactive() -> None:
    """Arrow-key selector to update any single credential in ~/.isaura/.env."""
    existing = _read_env_file(ENV_PATH)
    if not existing:
        console.print(
            "[red]No configuration file found.[/] "
            "Run any isaura command first to create it, then try again."
        )
        sys.exit(1)

    choices = [
        questionary.Choice(
            title=f"{key}  (current: {_mask(val) if _is_secret(key) else val})",
            value=key,
        )
        for key, val in existing.items()
    ]

    selected_key = questionary.select(
        "Which credential do you want to update?",
        choices=choices,
    ).ask()

    if selected_key is None:
        console.print("[dim]No changes made.[/]")
        return

    if _is_secret(selected_key):
        new_val = questionary.password(f"New value for {selected_key}:").ask()
    else:
        new_val = questionary.text(
            f"New value for {selected_key}:",
            default=existing[selected_key],
        ).ask()

    if new_val is None:
        console.print("[dim]No changes made.[/]")
        return

    confirmed = questionary.confirm(
        f"Are you sure you want to update {selected_key}?",
        default=False,
    ).ask()

    if confirmed:
        existing[selected_key] = new_val
        _write_env_file(ENV_PATH, existing)
        console.print(f"[green]{selected_key} updated successfully.[/]")
    else:
        console.print("[dim]No changes made.[/]")
