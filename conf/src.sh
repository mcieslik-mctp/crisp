#!/usr/bin/env bash
## workaround (silent if fails)
sudo su -c 'echo "kernel.shmmax = 31000000000" >> /etc/sysctl.conf; echo "kernel.shmall = 31000000000" >> /etc/sysctl.conf; /sbin/sysctl -p' > /dev/null 2>&1
## update path
export PATH=/code/bin:$PATH
