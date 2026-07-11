# AquaEsterNet-StructMine Online Demonstration Guide

## Start Now

Double-click once:

```text
install_online_ui.bat
```

The installer verifies the local UI, starts an outbound Cloudflare HTTPS tunnel without opening a router port, enables a random strong application password, and opens both the online page and `.cache/remote_ui/SHARE_WITH_EXPERTS.txt`.

Share the **Online URL** and **Access password** from that file only with invited experts. They can use the Run, AF Prediction, Candidates, Structure, History, and Doctor views from another computer in real time.

## Controls

| Double-click file | Purpose |
|---|---|
| `start_online_ui.bat` | Start or recover online access and open the current URL/share file |
| `online_ui_status.bat` | Inspect local health, tunnel health, processes, and current URL |
| `reset_online_password.bat` | Rotate the shared password and rebuild the online tunnel |
| `stop_online_ui.bat` | Close public access immediately while preserving the local UI and results |
| `uninstall_online_ui.bat` | Remove online auto-start and local access secrets without deleting the project |

## Demonstration Tunnel Limits

The default launcher uses a Cloudflare Quick Tunnel. This is suitable for a temporary expert demonstration, not an uptime-guaranteed hosted service:

- The URL may change after a tunnel restart; the share file is updated automatically.
- Computation remains on the owner's workstation, which must stay powered on, signed in, and online.
- Sleep, network loss, or unavailable WSL/GPU components interrupt corresponding operations.
- The shared password grants permission to start compute jobs, cancel tasks, and view existing results. Stop access or rotate the password after a demonstration.

Cloudflare also defines Quick Tunnels as testing and development facilities. A durable fixed address requires a named Tunnel, a domain, and preferably a Cloudflare Access allowlist:

- [Official Quick Tunnel documentation](https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/do-more-with-tunnels/trycloudflare/)
- [Official named Tunnel documentation](https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/get-started/create-remote-tunnel/)

## Security Boundary

- Streamlit continues to bind only to `127.0.0.1:8501`; no inbound firewall or router port is opened.
- External transport is HTTPS. Passwords are checked using a random salt and 600,000 PBKDF2-SHA256 iterations.
- Plaintext share details stay in the Git-ignored `.cache/remote_ui/` directory with a current-user ACL.
- Online mode disables server-side “open results folder” actions.
- Third-party public-MSA transfer still requires per-job confirmation; remote access does not change that privacy boundary.
- Structural predictions and WASE-screen rankings remain research triage and require experimental validation.

## Fixed Hostname Mode

Create a named Tunnel in Cloudflare, route a public hostname to `http://127.0.0.1:8501`, and obtain its Tunnel token. Then save the token for the current Windows user and install Token mode:

```powershell
[Environment]::SetEnvironmentVariable("AQUA_CLOUDFLARED_TOKEN", "<Tunnel token>", "User")

powershell -NoProfile -ExecutionPolicy Bypass -File scripts/windows/remote_ui_service.ps1 `
  -Action Install -Mode Token -PublicUrl https://your-host.example.com -OpenBrowser
```

The token is encrypted with current-user Windows DPAPI before it is stored in the ignored runtime directory. Configure a Cloudflare Access email allowlist as the first layer and retain the application password as a second layer.
