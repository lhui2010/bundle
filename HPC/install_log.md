

===

###  Install Cmake for hyphy

#### Motives

Upgrade Cmake

#### Problem met

```
Scanning dependencies of target CMakeLib
/bin/ld: warning: libcrypto.so.1.1, needed by /usr/lib64/libssl.so, may conflict with libcrypto.so.10
/bin/ld: lib/libcmcurl.a(openssl.c.o): undefined reference to symbol 'DSA_get0_key@@OPENSSL_1_1_0'
/usr/lib64/libcrypto.so.1.1: error adding symbols: DSO missing from command line
collect2: 错误：ld 返回 1
gmake[2]: *** [Utilities/cmcurl/curltest] 错误 1
gmake[1]: *** [Utilities/cmcurl/CMakeFiles/curltest.dir/all] 错误 2
gmake[1]: *** 正在等待未完成的任务....
```

#### Fix

```
cd /usr/lib64/
mv libcrypto.so.10.bak /usr/lib64/libcrypto.so.1.0.2k backup/
ldconfig -v |less
mv libcrypto.so backup/
ln -s libcrypto.so.1.1 libcrypto.so
```

2020年 10月 08日 星期四 11:27:55 CST

===
