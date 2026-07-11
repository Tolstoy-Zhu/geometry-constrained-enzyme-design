"""Compatibility entry point for the TM-align-capable environment doctor."""
from bur01_structmine.doctor import print_doctor, run_doctor


if __name__ == "__main__":
    print_doctor(run_doctor("configs/benchmark_tmtools.yaml", check_network=False))
