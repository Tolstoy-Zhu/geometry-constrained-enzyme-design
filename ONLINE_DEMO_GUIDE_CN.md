# AquaEsterNet-StructMine 在线演示说明

## 立即开始

首次双击：

```text
install_online_ui.bat
```

程序会完成四件事：

1. 确认本地网页 `http://127.0.0.1:8501` 正常运行。
2. 启动从本机到 Cloudflare 的出站 HTTPS 通道，不开放路由器端口。
3. 为在线工作台启用随机强密码。
4. 打开在线网页和 `.cache/remote_ui/SHARE_WITH_EXPERTS.txt`。

把分享文件中的 **Online URL** 和 **Access password** 发给受邀专家。专家可在其他电脑的现代浏览器中实时操作运行、AF预测、候选、结构、历史和诊断视图。

## 日常控制

| 双击文件 | 作用 |
|---|---|
| `start_online_ui.bat` | 启动或恢复在线访问，并打开最新网址与分享文件 |
| `online_ui_status.bat` | 查看本地服务、HTTPS 通道、进程和网址状态 |
| `reset_online_password.bat` | 更换共享密码并重新建立在线通道 |
| `stop_online_ui.bat` | 立即关闭公网访问，保留本地网页和全部结果 |
| `uninstall_online_ui.bat` | 移除在线自启动任务和本地访问密钥，不删除项目 |

## 当前演示通道的限制

默认入口使用 Cloudflare Quick Tunnel。它适合临时展示，但不是有可用性承诺的正式托管服务：

- 通道重启后网址可能变化；最新网址会自动写入分享文件。
- 所有计算仍由所有者电脑完成，因此电脑必须开机、登录 Windows 并联网。
- Windows 休眠、断网或 WSL/GPU 环境不可用时，专家无法完成相应任务。
- 共享密码允许受邀者启动计算、取消任务和查看已有结果。演示结束后应停止在线访问或重置密码。

Cloudflare 官方也将 Quick Tunnel 定义为测试和开发用途。长期使用固定网址时，应在 Cloudflare 账户中创建命名 Tunnel、绑定自己的域名并配置 Access 允许名单：

- [Quick Tunnel 官方说明](https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/do-more-with-tunnels/trycloudflare/)
- [命名 Tunnel 官方说明](https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/get-started/create-remote-tunnel/)

## 安全边界

- Streamlit 仍只监听 `127.0.0.1:8501`，Windows 防火墙和路由器不需要开放入站端口。
- 外部访问通过 HTTPS；访问密码使用随机盐和 `600,000` 次 PBKDF2-SHA256 派生后验证。
- 明文分享密码仅保存在被 Git 忽略的 `.cache/remote_ui/` 中，并限制为当前 Windows 用户访问。
- 在线模式禁用“打开结果目录”，防止远程点击在所有者电脑上弹出资源管理器。
- AlphaFold 公共 MSA 的第三方序列传输仍要求逐次确认；在线化没有改变该隐私边界。
- 本项目不是医疗或生产控制系统。专家操作产生的候选和评分仍需实验验证。

## 固定域名模式

固定域名需要先在 Cloudflare 控制台创建命名 Tunnel，将公开主机名指向 `http://127.0.0.1:8501`，并取得 Tunnel token。随后在当前 Windows 用户环境中保存 token，再安装 Token 模式：

```powershell
[Environment]::SetEnvironmentVariable("AQUA_CLOUDFLARED_TOKEN", "<Tunnel token>", "User")

powershell -NoProfile -ExecutionPolicy Bypass -File scripts/windows/remote_ui_service.ps1 `
  -Action Install -Mode Token -PublicUrl https://your-host.example.com -OpenBrowser
```

Token 会使用当前 Windows 用户的 DPAPI 加密后写入忽略目录，不进入源码。建议同时在 Cloudflare Access 中按专家邮箱建立允许名单；应用内共享密码保留为第二层保护。
