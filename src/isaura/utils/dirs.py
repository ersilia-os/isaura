"""Dir utils."""

import os
from pathlib import Path
from typing import Optional

from xdg import xdg_data_home


def get_workspace_path(override: Optional[Path] = None) -> Path:
    """Get path of the Olinda workspace.

    Args:
        override (Optional[Path], optional): _description_. Defaults to None.

    Returns:
        Path: Package root path
    """
    if override is None:
        workspace_path = xdg_data_home() / "eos"
    else:
        workspace_path = Path(override)
    # Create dir if not already present
    os.makedirs(workspace_path, exist_ok=True)
    return workspace_path
