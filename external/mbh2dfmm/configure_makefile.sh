#!/bin/sh

host=linux-gfortran
prefix=/usr

debugsym=false
openmpsym=true

for arg in "$@"; do
    case "$arg" in
    --prefix=*)
        prefix=`echo $arg | sed 's/--prefix=//'`
        ;;

    --host=*)
        host=`echo $arg | sed 's/--host=//'`
        ;;

    --enable-debug)
        debugsym=true;;
    --disable-debug)
        debugsym=false;;

    --enable-openmp)
        openmpsym=true;;
    --disable-openmp)
        openmpsym=false;;

    --list-systems)
	echo 'valid host names: '
	echo ' --host=macosx'
	echo ' --host=linux-gfortran'
	exit 0
	;;

    --help)
        echo 'usage: ./configure [options]'
        echo 'options:'
        echo '  --prefix=<path>: installation prefix'
        echo '  --host=<system-name>:  prefix'
	echo '  --list-systems: print valid host names'
        echo '  --enable-debug: include debug symbols'
        echo '  --disable-debug: do not include debug symbols'
	echo '  --enable-openmp: include openmp symbols'
	echo '  --disable-openmp: do not include openmp symbols'
        echo 'all invalid options are silently ignored'
        exit 0
        ;;
    esac
done

echo 'generating makefile ...'
echo "PREFIX = $prefix" >Makefile
echo "HOST = $host" >>Makefile
if $debugsym; then
    echo 'DBGYN = yes' >>Makefile
else
    echo 'DBGYN = no' >>Makefile
fi
if $openmpsym; then
    echo 'OPENMPYN = yes' >>Makefile
else
    echo 'OPENMPYN = no' >>Makefile
fi
cat Makefile.in >>Makefile
echo 'configuration complete, type make to build.'
