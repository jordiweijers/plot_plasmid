import logging
import os
import subprocess
from typing import List

def setup_logging(log_file, name="main_logger", console=True):
    """ Setup up logging."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    if logger.hasHandlers():
        logger.handlers.clear()
	# Console handler
    if console:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter("%(message)s")
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
    # File handler
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)-8s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
	)
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)
    return logger

def run_command(cmd, log_path, cpus=1):
    log_dir = os.path.dirname(log_path)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": str(cpus),
        "OPENBLAS_NUM_THREADS": str(cpus),
        "MKL_NUM_THREADS": str(cpus),
        "NUMEXPR_NUM_THREADS": str(cpus),
        "VECLIB_MAXIMUM_THREADS": str(cpus),
        "BLIS_NUM_THREADS": str(cpus),
    })
    with open(log_path, "w") as log:
        process = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            env=env
        )
    if process.returncode != 0:
        raise RuntimeError(f"External command failed (see {log_path})")