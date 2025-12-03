# 使用方法：
#   source ~/vpn.sh on [host:port]   # 开启代理（默认 127.0.0.1:7897）
#   source ~/vpn.sh off              # 关闭代理
#   source ~/vpn.sh status           # 查看当前 _proxy 变量

case "$1" in
  on)
    HOSTPORT="${2:-127.0.0.1:7898}"
    export http_proxy="http://$HOSTPORT"
    export https_proxy="http://$HOSTPORT"
    export all_proxy="socks5://$HOSTPORT"
    # 同时设置大写，兼容某些程序
    export HTTP_PROXY="$http_proxy"
    export HTTPS_PROXY="$https_proxy"
    export ALL_PROXY="$all_proxy"
    # 建议直连的地址（按需增删内网网段）
    export no_proxy="localhost,127.0.0.1,::1"
    export NO_PROXY="$no_proxy"
    echo "[vpn] Proxy ON -> $HOSTPORT"
    echo "[vpn] Test: curl -I https://google.com  或  curl -s https://ipinfo.io/ip"
    ;;
  off)
    unset http_proxy https_proxy all_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY
    unset no_proxy NO_PROXY
    echo "[vpn] Proxy OFF"
    ;;
  status)
    env | grep -i _proxy || echo "[vpn] no *_proxy set"
    ;;
  *)
    echo "Usage: source ~/vpn.sh {on|off|status} [host:port]"
    ;;
esac
