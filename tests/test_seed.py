# tests/test_seed.py
import os
import random

import numpy as np
import torch
from seed import seed_everything


def test_seed_everything():
    seed_value = 42

    # 1. Set the seed
    seed_everything(seed_value)

    # 2. Generate random numbers
    py_rand1 = random.random()
    np_rand1 = np.random.rand()
    torch_rand1 = torch.rand(1).item()

    # 3. Reset the seed to the SAME value
    seed_everything(seed_value)

    # 4. Generate random numbers again
    py_rand2 = random.random()
    np_rand2 = np.random.rand()
    torch_rand2 = torch.rand(1).item()

    # 5. Assertions: They must be exactly equal
    assert py_rand1 == py_rand2, "Python random module not seeded correctly"
    assert np_rand1 == np_rand2, "Numpy not seeded correctly"
    assert torch_rand1 == torch_rand2, "PyTorch not seeded correctly"
    assert os.environ["PYTHONHASHSEED"] == str(seed_value), "Env var not set"
