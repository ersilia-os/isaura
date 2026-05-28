import pyarrow as pa
import pyarrow.csv as pa_csv
from rich.console import Console
from rich.progress import (
  Progress, ProgressColumn, SpinnerColumn, TextColumn,
  BarColumn, MofNCompleteColumn, TaskProgressColumn,
  TimeRemainingColumn, TimeElapsedColumn,
)
from rich.progress_bar import ProgressBar
from rich.table import Table

from isaura.utils import rss_mb


class StreamingCsvSink:
  """Context manager that writes DataFrames to a CSV file incrementally.

  Writes the header only once (on the first write_table call) and appends
  subsequent batches without re-writing it. Used to stream large result sets
  to disk without holding everything in memory.

  Args:
      path: File path to write the CSV to.
  """

  def __init__(self, path):
    self._path = path
    self._fp = None
    self._header_written = False
    self.rows_written = 0

  def __enter__(self):
    self._fp = open(self._path, "wb")
    return self

  def __exit__(self, et, ev, tb):
    self.close()

  def close(self):
    if self._fp is not None and not self._fp.closed:
      self._fp.close()

  def write_table(self, table):
    """Append a PyArrow Table to the CSV, writing the header only on the first call."""
    if table is None or table.num_rows == 0:
      return
    opts = pa_csv.WriteOptions(include_header=not self._header_written, quoting_style="none", quoting_header="none")
    pa_csv.write_csv(table, self._fp, write_options=opts)
    self.rows_written += table.num_rows
    self._header_written = True

  def write_batch(self, batch):
    """Append a PyArrow RecordBatch to the CSV."""
    self.write_table(pa.Table.from_batches([batch]))

  def write_df(self, df):
    """Append a pandas DataFrame to the CSV."""
    if df is None or df.empty:
      return
    self.write_table(pa.Table.from_pandas(df, preserve_index=False))


class _PulseBarColumn(ProgressColumn):
  """Rich progress column that renders a pulsing bar while a task is in progress.

  Switches to a solid filled bar when the task is complete.
  """

  def __init__(
    self, bar_width=40,
    style="rgb(36,26,58)",
    complete_style="bold rgb(56,189,248)",
    finished_style="bold rgb(34,197,94)",
    pulse_style="bold rgb(244,114,182)",
  ):
    super().__init__()
    self.bar_width = int(bar_width)
    self.style = style
    self.complete_style = complete_style
    self.finished_style = finished_style
    self.pulse_style = pulse_style

  def render(self, task):
    return ProgressBar(
      total=task.total, completed=task.completed,
      width=max(1, self.bar_width), pulse=not task.finished,
      animation_time=task.get_time(),
      style=self.style, complete_style=self.complete_style,
      finished_style=self.finished_style, pulse_style=self.pulse_style,
    )


class ReadProgress:
  """Context manager that displays a live Rich progress bar during a read operation.

  Shows files scanned, rows emitted, pending inputs, memory usage, and
  elapsed/remaining time. Designed to be used with the `with` statement and
  updated incrementally via update().

  Args:
      total_inputs: Total number of molecule inputs being queried (used to size the bar).
      console: Rich Console to render to. Defaults to a new stderr console.
  """

  def __init__(self, total_inputs=None, console=None):
    self.total_inputs = total_inputs
    self.console = console or Console(force_terminal=True, stderr=True)
    self.progress = None
    self.task_id = None
    self.files_done = 0
    self.files_total = 0
    self.found_rows = 0
    self.emitted_rows = 0
    self.unresolved = total_inputs or 0
    self._tick = 0
    self._phrases = ["scanning", "matching", "streaming", "finalizing"]

  def __enter__(self):
    self.progress = Progress(
      SpinnerColumn(spinner_name="aesthetic"),
      TextColumn("[bold magenta]{task.fields[phase]}[/]"),
      _PulseBarColumn(bar_width=28, style="rgb(40,28,58)",
        complete_style="bold rgb(56,189,248)",
        finished_style="bold rgb(34,197,94)",
        pulse_style="bold rgb(244,114,182)"),
      TaskProgressColumn(text_format="[bold bright_cyan]{task.percentage:>3.0f}%[/]"),
      MofNCompleteColumn(),
      TextColumn("[cyan]files[/] {task.fields[files_done]}/{task.fields[files_total]}"),
      TextColumn("[green]emitted[/] {task.fields[emitted_rows]}"),
      TextColumn("[yellow]pending[/] {task.fields[unresolved]}"),
      TextColumn("[magenta]rss[/] {task.fields[rss]}MB"),
      TimeElapsedColumn(), TimeRemainingColumn(),
      console=self.console, transient=False,
      redirect_stdout=False, redirect_stderr=False,
      refresh_per_second=10, expand=True,
    )
    self.progress.__enter__()
    total = self.total_inputs if self.total_inputs is not None else 100
    self.task_id = self.progress.add_task(
      self._render("starting"), total=total, completed=0,
      phase=self._render_phase("starting"), files_done=0, files_total=0,
      emitted_rows=0, unresolved=self.unresolved, rss=f"{rss_mb():.0f}",
    )
    return self

  def __exit__(self, et, ev, tb):
    if self.progress is not None:
      self.progress.__exit__(et, ev, tb)

  def _render(self, stage):
    phrase = self._phrases[self._tick % len(self._phrases)]
    return f"{phrase} {stage}"

  def _render_phase(self, stage):
    return self._render(stage).replace("_", " ")

  def update(self, stage=None, files_done=None, files_total=None,
             found_rows=None, emitted_rows=None, unresolved=None):
    """Update the progress bar with the latest scan state."""
    if files_done is not None:
      self.files_done = files_done
    if files_total is not None:
      self.files_total = files_total
    if found_rows is not None:
      self.found_rows = found_rows
    if emitted_rows is not None:
      self.emitted_rows = emitted_rows
    if unresolved is not None:
      self.unresolved = unresolved
    self._tick += 1
    if self.progress is not None and self.task_id is not None:
      completed = self.found_rows
      total = self.total_inputs if self.total_inputs is not None else max(100, completed or 0)
      if completed > total:
        total = completed
      self.progress.update(
        self.task_id, description=self._render(stage or "working"),
        total=total, completed=completed,
        phase=self._render_phase(stage or "working"),
        files_done=self.files_done, files_total=self.files_total,
        emitted_rows=self.emitted_rows, unresolved=self.unresolved,
        rss=f"{rss_mb():.0f}",
      )


def track_write_progress(rows, total=None, description="Writing rows", console=None):
  """Wrap an iterable of rows with a Rich progress bar and yield each row through.

  Args:
      rows: Iterable of rows to pass through.
      total: Total row count for a determinate bar. Uses a spinner if None.
      description: Label shown next to the progress bar.
      console: Rich Console to render to.

  Yields:
      Each row from the input iterable unchanged.
  """
  if total is None:
    progress = Progress(
      SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
      TextColumn("{task.completed} rows"), TimeElapsedColumn(),
      console=console, transient=True,
    )
  else:
    progress = Progress(
      SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
      BarColumn(), TextColumn("{task.completed}/{task.total}"),
      TimeElapsedColumn(), TimeRemainingColumn(),
      console=console, transient=True,
    )
  with progress:
    task_id = progress.add_task(description, total=total)
    for row in rows:
      yield row
      progress.advance(task_id, 1)


def spinner(message, fn, *args, **kwargs):
  """Run fn(*args, **kwargs) while displaying a Rich spinner with message. Returns the result."""
  c = Console()
  with c.status(message, spinner="dots"):
    result = fn(*args, **kwargs)
  return result


def make_table(title, cols, rows):
  """Build a Rich Table from a list of column specs and row dicts.

  Args:
      title: Table title string.
      cols: List of dicts with "name", optional "justify", and optional "style" keys.
      rows: List of dicts mapping column "key" values to cell content.

  Returns:
      A Rich Table ready to print.
  """
  t = Table(title=title)
  for c in cols:
    t.add_column(c["name"], justify=c.get("justify", "left"), style=c.get("style", ""))
  for r in rows:
    t.add_row(*[str(r.get(c["key"], "")) for c in cols])
  return t
