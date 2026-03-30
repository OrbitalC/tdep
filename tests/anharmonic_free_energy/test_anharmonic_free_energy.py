import numpy as np
import pytest
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"
cases = [
    "Silicon_1600K_classical",
    "Neon_14K_quantum",
]
output_file = "outfile.anharmonic_thermodynamics"
# Seems to be some noise due to ordering of sums.
rtol = 1e-6
atol = 1e-8


def _read_file(file):
    values = []
    with open(file, encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            row = np.fromstring(stripped, sep=" ")
            if row.size:
                values.extend(row.tolist())
    return np.asarray(values, dtype=float)


@pytest.mark.parametrize("case", cases)
def test_output(case):
    file_ref = folder / f"{case}_reference"
    file_new = parent / case / output_file
    data_ref = _read_file(file_ref)
    data_new = _read_file(file_new)

    assert data_ref.shape == data_new.shape, file_new.absolute()
    np.testing.assert_allclose(
        data_ref,
        data_new,
        rtol=rtol,
        atol=atol,
        err_msg=str(file_new.absolute()),
    )


if __name__ == "__main__":
    for case in cases:
        test_output(case)
