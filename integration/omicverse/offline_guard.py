"""Test-only network guard: block HTTP, DNS and INET socket traffic.

Does not replace OmicVerse modules, numerical functions or dependency checks.
No network settings are changed outside this process. Unix-domain sockets and
loopback are allowed only when allow_loopback=True (notebook kernel transport).
"""
import ipaddress
import socket


def install(*, allow_loopback=False):
    blocked = []
    def allowed(address):
        if not allow_loopback or not isinstance(address, tuple): return False
        try: return ipaddress.ip_address(address[0]).is_loopback
        except ValueError: return address[0] == 'localhost'
    def wrap(old, label):
        def guarded(self, address, *args, **kwargs):
            if self.family == socket.AF_UNIX or allowed(address):
                return old(self, address, *args, **kwargs)
            blocked.append(label)
            raise OSError('network disabled during Scissor integration validation')
        return guarded
    socket.socket.connect = wrap(socket.socket.connect, 'socket.connect')
    socket.socket.connect_ex = wrap(socket.socket.connect_ex, 'socket.connect_ex')
    old_dns = socket.getaddrinfo
    def dns(host, *args, **kwargs):
        if allow_loopback and (host in ('localhost', '127.0.0.1', '::1', None)):
            return old_dns(host, *args, **kwargs)
        blocked.append('DNS')
        raise OSError('DNS disabled during Scissor integration validation')
    socket.getaddrinfo = dns
    def sendto(self, *args, **kwargs):
        blocked.append('socket.sendto')
        raise OSError('UDP disabled during Scissor integration validation')
    socket.socket.sendto = sendto
    if hasattr(socket.socket, 'sendmsg'):
        socket.socket.sendmsg = sendto
    import requests
    def request(self, *args, **kwargs):
        blocked.append('requests')
        raise OSError('HTTP disabled during Scissor integration validation')
    requests.Session.request = request
    return blocked
