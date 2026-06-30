import pyarrow as pa
import pyarrow.csv as pa_csv
from rich.console import Console
from rich.progress import (
  Progress, SpinnerColumn, TextColumn,
  BarColumn, TaskProgressColumn,
  TimeRemainingColumn, TimeElapsedColumn,
)
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


class ReadProgress:
  """Context manager that displays a live Rich progress bar during a read operation.

  Shows files processed, row count, and elapsed/remaining time.
  Designed to be used with the `with` statement and updated incrementally via update().

  Args:
      total_inputs: Total number of molecule inputs being queried (used to size the bar).
      console: Rich Console to render to. Defaults to a new stderr console.
      description: Label shown next to the progress bar.
  """

  def __init__(self, total_inputs=None, console=None, description="Reading", writing_description=None):
    self.total_inputs = total_inputs
    self.console = console or Console(force_terminal=True, stderr=True)
    self.description = description
    # Shown once the read flips from scanning files to sorting/writing the
    # output. Falls back to the reading label so callers without a distinct
    # write phase (e.g. pull) keep a single description throughout.
    # A distinct write phase exists only when the caller passes a writing label
    # (i.e. a CSV is being written). Without it (pull / read without --output)
    # we keep a single morphing bar.
    self._has_write_phase = writing_description is not None
    self.writing_description = writing_description or description
    self.stage = None
    self.progress = None
    self.task_id = None
    self.write_task_id = None
    self.found_rows = 0
    self.files_done = 0
    self.files_total = 0
    self.emitted_rows = 0
    self.unresolved = total_inputs or 0

  def __enter__(self):
    total = self.total_inputs
    self.progress = Progress(
      SpinnerColumn(),
      TextColumn("[progress.description]{task.description}"),
      BarColumn(),
      TaskProgressColumn(text_format="[bold bright_cyan]{task.percentage:>3.0f}%[/]"),
      TextColumn("{task.completed}/{task.total}"),
      TextColumn("[cyan]files[/] {task.fields[files_done]}/{task.fields[files_total]}"),
      TimeElapsedColumn(),
      TimeRemainingColumn(),
      console=self.console, transient=False,
      redirect_stdout=False, redirect_stderr=False,
      refresh_per_second=10, expand=True,
    )
    self.progress.__enter__()
    self.task_id = self.progress.add_task(
      self.description, total=total, completed=0,
      files_done=0, files_total=0,
    )
    return self

  def __exit__(self, et, ev, tb):
    if self.progress is not None:
      self.progress.__exit__(et, ev, tb)

  def update(self, stage=None, files_done=None, files_total=None,
             found_rows=None, emitted_rows=None, unresolved=None):
    """Update the progress bar with the latest scan state."""
    if stage is not None:
      self.stage = stage
    if found_rows is not None:
      self.found_rows = found_rows
    if files_done is not None:
      self.files_done = files_done
    if files_total is not None:
      self.files_total = files_total
    if emitted_rows is not None:
      self.emitted_rows = emitted_rows
    if unresolved is not None:
      self.unresolved = unresolved
    if self.progress is None or self.task_id is None:
      return

    def _total_for(completed):
      t = self.total_inputs if self.total_inputs is not None else max(100, completed or 0)
      return max(t, completed)

    if self.stage == "emitting" and self._has_write_phase:
      # Reading is done. Freeze the read bar at 100% (kept on screen because the
      # Progress is non-transient) and drive a SEPARATE write bar below it, so
      # both the read time and the CSV-write time stay visible side by side.
      read_total = _total_for(self.found_rows)
      self.progress.update(
        self.task_id, completed=read_total, total=read_total,
        files_done=self.files_done, files_total=self.files_total,
      )
      if self.write_task_id is None:
        self.write_task_id = self.progress.add_task(
          self.writing_description, total=_total_for(self.emitted_rows), completed=0,
          files_done=self.files_done, files_total=self.files_total,
        )
      self.progress.update(
        self.write_task_id, completed=self.emitted_rows, total=_total_for(self.emitted_rows),
        files_done=self.files_done, files_total=self.files_total,
      )
    elif self.stage == "emitting":
      # No distinct write phase (pull / read without --output): single bar that
      # tracks emitted rows once the scan is done.
      self.progress.update(
        self.task_id, completed=self.emitted_rows, total=_total_for(self.emitted_rows),
        description=self.writing_description,
        files_done=self.files_done, files_total=self.files_total,
      )
    else:
      self.progress.update(
        self.task_id, completed=self.found_rows, total=_total_for(self.found_rows),
        description=self.description,
        files_done=self.files_done, files_total=self.files_total,
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
