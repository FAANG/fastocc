import os

import pytest
from click.testing import CliRunner
from fastocc.fastocc import fastocc


inputs = {
    "bed": "data/regions.bed",
    "bam": "data/reads.bam",
    "opts": "data/opts.json"
}


outputs = {
    "fragments": "data/fragments.h5wig",
    "nfrs": "data/nfrs.bed",
    "occ": "data/occ.bw",
    "bg": "data/fragments.bedgraph"
}


@pytest.fixture(autouse=True, scope="class")
def cleanup_files():
    yield
    for file in outputs.values():
        if os.path.exists(file):
            os.remove(file)


class TestCli:
    def test_fastocc_call(self):
        if os.path.exists("test.h5wig"):
            os.remove("test.h5wig")

        runner = CliRunner()
        result = runner.invoke(fastocc, ["call", "--bed", inputs["bed"], "--bam", inputs["bam"], "--output", outputs["fragments"]])

        assert os.path.exists(outputs["fragments"])
        assert result.exit_code == 0

    def test_fastocc_nfr(self):
        if not os.path.exists(outputs["fragments"]):
            pytest.skip("file 'test.h5wig' required to test `fastocc nfr`")

        runner = CliRunner()
        result = runner.invoke(fastocc, ["nfr", "--fragments", outputs["fragments"], "--output", outputs["nfrs"]])

        assert os.path.exists(outputs["nfrs"])
        assert result.exit_code == 0

    def test_fastocc_export(self):
        if not os.path.exists(outputs["fragments"]):
            pytest.skip("file 'test.h5wig' required to test `fastocc export`")

        runner = CliRunner()
        result = runner.invoke(fastocc, ["export", "--fragments", outputs["fragments"], "--output", outputs["occ"]])

        assert os.path.exists(outputs["occ"])
        assert result.exit_code == 0

    def test_fastocc_dump(self):
        if not os.path.exists(outputs["fragments"]):
            pytest.skip("file 'test.h5wig' required to test `fastocc export`")

        runner = CliRunner()
        result = runner.invoke(fastocc, ["dump", "--fragments", outputs["fragments"], "--bed", inputs["bed"], "--output", outputs["bg"]])

        assert os.path.exists(outputs["bg"])
        assert result.exit_code == 0

def test_fastocc_nps():
    if os.path.exists(outputs["fragments"]):
        os.remove(outputs["fragments"])

    runner = CliRunner()
    result = runner.invoke(fastocc, ["call", "--bed", inputs["bed"], "--bam", inputs["bam"], "--output", outputs["fragments"], "--nps"])

    assert os.path.exists(outputs["fragments"])
    assert result.exit_code == 0


def test_fastocc_opts():
    if os.path.exists(outputs["fragments"]):
        os.remove(outputs["fragments"])

    runner = CliRunner()
    result = runner.invoke(fastocc, ["call", "--bed", inputs["bed"], "--bam", inputs["bam"], "--output", outputs["fragments"], "--options", inputs["opts"]])

    assert os.path.exists(outputs["fragments"])
    assert result.exit_code == 0
