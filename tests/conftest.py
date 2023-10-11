import os
import pytest


@pytest.fixture(autouse=True)
def run_around_tests(capsys):
    yield

    with capsys.disabled():
        os.system('free -h')
        os.system('df -h')
