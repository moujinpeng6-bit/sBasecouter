"""An OS-released lock for output directories with resumable state."""

from contextlib import contextmanager
import os


@contextmanager
def output_lock(out):
    with (out / ".writer.lock").open("a+b") as lock:
        try:
            if os.name == "nt":
                import msvcrt
                lock.seek(0)
                msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise ValueError(f"Another task is using the output directory: {out}") from exc
        yield
