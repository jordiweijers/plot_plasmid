import logging
import os
import subprocess
from tqdm import tqdm
from typing import List, Callable, Iterable
import multiprocessing

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
    
def run_in_parallel(function: Callable, args_list: Iterable[List], cpus: int) -> List:
	"""
	Run a given function in parallel over available CPUs.
	Args:
		function (Callable): The function to call.
		args_list(Iterable[List]): Iterabale of argument lists for each function call.
		cpus (int): Number ofCPUs to use.
	Returns:
		List: List of return values from each function call.
	"""
	if cpus <= 1:
		results = []
		for args in tqdm(args_list, desc="Processing"):
			results.append(function(*args))
		return results
	for item in args_list:
		try:
			iter(item)
		except TypeError as e:
			raise TypeError(f"Each element of args_list must be an iterable of arguments. Got: {item}")
	manager = multiprocessing.Manager()
	results = manager.list()
	pbar = tqdm(total=len(args_list), desc="Processing")
	def callback(result):
		results.append(result)
		pbar.update()
	with multiprocessing.Pool(cpus) as pool:
		for args in args_list:
			pool.apply_async(function, args=args, callback=callback)
		pool.close()
		pool.join()
	pbar.close()
	return list(results)