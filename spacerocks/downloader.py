import os.path
import sys
from concurrent.futures import ThreadPoolExecutor
import signal
from functools import partial
from threading import Event
from typing import Iterable
from urllib.request import urlopen
import contextlib

from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TaskID,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

done_event = Event()

def handle_sigint(signum, frame):
    done_event.set()

signal.signal(signal.SIGINT, handle_sigint)


def download(urls: Iterable[str], dest_dir: str):
    """Download multuple files to the given directory."""

    progress = Progress(TextColumn("[bold blue]{task.fields[filename]}", justify="right"),
                                   BarColumn(bar_width=None),
                                   "[progress.percentage]{task.percentage:>3.1f}%",
                                   "•",
                                   DownloadColumn())

    def copy_url(task_id: TaskID, url: str, path: str) -> None:
        """Copy data from a url to a local file."""
        with contextlib.closing(urlopen(url)) as response:
            # This will break if the response doesn't contain content length
            progress.update(task_id, total=int(response.info()["Content-length"]))
            with open(path, "wb") as dest_file:
                progress.start_task(task_id)
                for data in iter(partial(response.read, 32768), b""):
                    dest_file.write(data)
                    progress.update(task_id, advance=len(data))
                    if done_event.is_set():
                        return

    with progress:
        dummy = []
        with ThreadPoolExecutor(max_workers=4) as pool:
            for url in urls:
                filename = url.split("/")[-1]
                dest_path = os.path.join(dest_dir, filename)
                task_id = progress.add_task("download", filename=filename, start=False)
                dummy.append(pool.submit(copy_url, task_id, url, dest_path))

        _ = [d.result() for d in dummy]