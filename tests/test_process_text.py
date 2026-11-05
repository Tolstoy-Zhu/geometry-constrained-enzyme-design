import codecs
import subprocess
import sys

import pytest

from bur01_structmine.process_text import capture_process_text, decode_process_output, read_process_log


MESSAGE = "WSL2 无法启动，因为此计算机上未启用虚拟化。\r\nHCS_E_HYPERV_NOT_INSTALLED\r\n"


@pytest.mark.parametrize(
    ("prefix", "encoding"),
    [
        (b"", "utf-16-le"),
        (codecs.BOM_UTF16_LE, "utf-16-le"),
        (codecs.BOM_UTF16_BE, "utf-16-be"),
        (b"", "utf-8"),
        (codecs.BOM_UTF8, "utf-8"),
    ],
)
def test_wsl_and_linux_output_keep_the_original_text(prefix, encoding):
    assert decode_process_output(prefix + MESSAGE.encode(encoding)) == MESSAGE


def test_probe_keeps_exit_status_and_decodes_streams_independently(monkeypatch):
    def failed_probe(command, *, env, **kwargs):
        assert env["WSL_UTF8"] == "1" and env["PYTHONIOENCODING"] == "utf-8"
        return subprocess.CompletedProcess(command, 7, "工具输出".encode(), MESSAGE.encode("utf-16-le"))

    monkeypatch.setattr(subprocess, "run", failed_probe)
    with pytest.raises(subprocess.CalledProcessError) as failure:
        capture_process_text(["wsl.exe"], check=True)
    assert failure.value.returncode == 7
    assert failure.value.stdout == "工具输出"
    assert failure.value.stderr == MESSAGE


def test_redirected_python_job_prints_utf8_even_on_windows(monkeypatch):
    monkeypatch.setenv("PYTHONIOENCODING", "ascii")
    command = [sys.executable, "-c", "print('候选序列')"]
    assert capture_process_text(command, check=True).stdout.strip() == "候选序列"


def test_existing_log_with_utf8_commands_and_utf16_wsl_diagnostics(tmp_path):
    prefix, suffix = "开始结构挖掘\r\n$ wsl.exe\r\n", "[CMD] 第二步\r\n"
    path = tmp_path / "ui_process.log"
    data = prefix.encode() + MESSAGE.encode("utf-16-le") + suffix.encode()
    path.write_bytes(data)
    assert read_process_log(path) == prefix + MESSAGE + suffix
    assert path.read_bytes() == data


def test_workbench_shows_same_wsl_cause_in_run_af_and_diagnostics(monkeypatch, tmp_path):
    from streamlit.testing.v1 import AppTest

    from bur01_structmine.config import project_root
    from bur01_structmine.doctor import DoctorCheck, DoctorReport

    report = DoctorReport(
        "test",
        "test",
        [
            DoctorCheck(
                "wsl_distribution",
                "wsl",
                "pass",
                "Ubuntu-Foldseek: installed",
                {"distribution": "Ubuntu-Foldseek"},
            ),
            DoctorCheck("wsl_gpu", "prediction", "fail", MESSAGE, {}),
            DoctorCheck("colabfold", "prediction", "warn", "ColabFoldUnavailable: " + MESSAGE, {}),
        ],
    )
    monkeypatch.delenv("AQUA_TEST_MODE", raising=False)
    monkeypatch.setenv("AQUA_REMOTE_ACCESS_PATH", str(tmp_path / "access.json"))
    monkeypatch.setenv("AQUA_RUN_ROOT", str(tmp_path / "runs"))
    monkeypatch.setenv("AQUA_PREDICTION_ROOT", str(tmp_path / "predictions"))
    monkeypatch.setenv("AQUA_NR_RUN_ROOT", str(tmp_path / "nr_runs"))
    monkeypatch.setattr("bur01_structmine.doctor.run_doctor", lambda *args, **kwargs: report)

    app = AppTest.from_file(str(project_root() / "ui/app.py"), default_timeout=60).run()
    assert not app.exception
    tabs = {tab.label: tab for tab in app.tabs}
    run_message = next(
        error.value for error in tabs["运行"].error if "HCS_E_HYPERV_NOT_INSTALLED" in error.value
    )
    assert "虚拟化" in run_message
    assert [warn.value for warn in tabs["AF预测"].warning].count(run_message) == 1
    diagnostics = tabs["诊断"].dataframe[0].value.set_index("check")
    assert diagnostics.loc["wsl_gpu", "message"] == run_message
    assert "尚未安装" not in diagnostics.loc["wsl_distribution", "message"]
