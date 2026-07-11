from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_persistent_ui_controller_has_required_safety_and_recovery_controls() -> None:
    controller = (ROOT / "scripts" / "windows" / "ui_service.ps1").read_text(encoding="utf-8-sig")

    assert '"http://127.0.0.1:$Port"' in controller
    assert '/_stcore/health' in controller
    assert '-AtLogOn' in controller
    assert '-RepetitionInterval (New-TimeSpan -Minutes 1)' in controller
    assert '-RestartCount 999' in controller
    assert '--server.fileWatcherType=none' in controller
    assert 'Test-ProjectUiProcess' in controller
    assert 'AQUA_UI_SERVICE' in controller


def test_windows_ui_entry_points_use_the_persistent_controller() -> None:
    expected_actions = {
        "start_ui.bat": "-Action Start",
        "install_ui_autostart.bat": "-Action Install",
        "stop_ui.bat": "-Action Stop",
        "ui_status.bat": "-Action Status",
        "uninstall_ui_autostart.bat": "-Action Uninstall",
    }

    for filename, action in expected_actions.items():
        content = (ROOT / filename).read_text(encoding="utf-8-sig")
        assert "scripts\\windows\\ui_service.ps1" in content
        assert action in content


def test_web_workbench_reports_persistent_service_mode() -> None:
    app = (ROOT / "ui" / "app.py").read_text(encoding="utf-8")

    assert 'os.environ.get("AQUA_UI_SERVICE") == "1"' in app
    assert "后台常驻 · 自动恢复" in app
    assert "Persistent local service · auto-recovery" in app


def test_online_ui_controller_keeps_origin_local_and_requires_authentication() -> None:
    controller = (ROOT / "scripts" / "windows" / "remote_ui_service.ps1").read_text(encoding="utf-8-sig")

    assert '$OriginUrl = "http://127.0.0.1:8501"' in controller
    assert "Rfc2898DeriveBytes" in controller
    assert "600000" in controller
    assert "trycloudflare" in controller
    assert "AQUA_CLOUDFLARED_TOKEN" in controller
    assert "Protect-SecretFile" in controller
    assert '-RepetitionInterval (New-TimeSpan -Minutes 1)' in controller


def test_online_ui_entry_points_use_the_remote_controller() -> None:
    expected_actions = {
        "install_online_ui.bat": "-Action Install",
        "start_online_ui.bat": "-Action Start",
        "stop_online_ui.bat": "-Action Stop",
        "online_ui_status.bat": "-Action Status",
        "reset_online_password.bat": "-Action Restart",
        "uninstall_online_ui.bat": "-Action Uninstall",
    }

    for filename, action in expected_actions.items():
        content = (ROOT / filename).read_text(encoding="utf-8-sig")
        assert "scripts\\windows\\remote_ui_service.ps1" in content
        assert action in content
