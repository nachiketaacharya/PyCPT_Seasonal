from __future__ import print_function
import platform
import subprocess as sp
import sys, os

if platform.system() == "Windows":
    had_to_get_mingw = 0
    had_to_get_sevenzip = 0
    try:
        print('Checking for 7zip... ', end='')
        sys.stdout.flush()
        f = open('7z_install_output.txt', 'w')
        sp.call(['7z'], stdout=f)
        f.close()
        print('     finished')
        had_to_get_sevenzip = 0
    except:
        had_to_get_sevenzip = 1
        print("Unexpected error:", sys.exc_info()[0])
        try:
            print('Installing 7zip... ', end='')
            sys.stdout.flush()
            sp.call(['curl', '-s', 'https://www.7-zip.org/a/7z1900.exe', '--output', './7z1900.exe'])
            sp.call(['7z1900.exe', '/S', '/D="C:\\Program Files\\7-Zip"'])
            path = os.environ.get('PATH')
            sp.call(['setx', 'PATH', '{};C:\\Program Files\\7-Zip'.format(path), '/M'])
            print('     finished')
            print('Dont forget to restart your shell so path takes effect!')
            try:
                print('Checking for 7zip again... ', end='')
                sys.stdout.flush()
                f = open('7z_install_output.txt', 'w')
                sp.call(['7z'],env=dict(os.environ, PATH='{};C:\\Program Files\\7-Zip'.format(path)), stdout=f)
                f.close()
                print('     finished')
            except:
                print("Unexpected error, go download 7zip by hand:", sys.exc_info()[0])
                sys.exit()
        except:
            print("Unexpected error, maybe run as administrator?:", sys.exc_info()[0])
            sys.exit()



    try:
        print('Checking for MinGW32-make... ', end='')
        f = open('check.txt', 'w')
        sp.call(['mingw32-make'], stdout=f, stderr=f)
        f.close()
        print('    finished')
        had_to_get_mingw = 0
    except:
        print("couldnt find:", sys.exc_info()[0])
        had_to_get_mingw= 1
        try:
            f = open('garbage2.txt','w')
            sp.call(['del', '/s', '/q', 'MinGW'], shell=True, stdout=f)
            sp.call(['rmdir', '/s', '/q', 'MinGW'], shell=True)
            sp.call(['del', '/s', '/q', 'mingw_install_output.txt'], shell=True, stdout=f)
            sp.call(['mkdir', 'MinGW'], shell=True, stdout=f)
            f.close()
            progress = 0
            print('fetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://newcontinuum.dl.sourceforge.net/project/mingw/MinGW/Base/binutils/binutils-2.28/binutils-2.28-1-mingw32-bin.tar.xz", '--output', './MinGW/binutils-2.28-1-mingw32-bin.tar.xz'])
            print('\rFetching MinGW...  {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', 'https://cfhcable.dl.sourceforge.net/project/mingw/MinGW/Base/mingwrt/mingwrt-5.0.2/mingwrt-5.0.2-mingw32-dev.tar.xz', '--output', './MinGW/mingwrt-5.0.2-dev.tar.xz'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', 'https://phoenixnap.dl.sourceforge.net/project/mingw/MinGW/Base/mingwrt/mingwrt-5.0.2/mingwrt-5.0.2-mingw32-dll.tar.xz', '--output', './MinGW/mingwrt-5.0.2-dll.tar.xz'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://cfhcable.dl.sourceforge.net/project/mingw/MinGW/Base/w32api/w32api-3.17/w32api-3.17-2-mingw32-dev.tar.lzma", '--output', './MinGW/w32api-3.17-2-mingw32-dev.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', 'https://pilotfiber.dl.sourceforge.net/project/mingw/MinGW/Base/mpc/mpc-1.0.1-2/mpc-1.0.1-2-mingw32-dll.tar.lzma', '--output', './MinGW/mpc-1.0.1-2-mingw32-dll.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', 'https://phoenixnap.dl.sourceforge.net/project/mingw/MinGW/Base/mpfr/mpfr-3.1.2-2/mpfr-3.1.2-2-mingw32-dll.tar.lzma', '--output', './MinGW/mpfr-3.1.2-2-mingw32-dll.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', 'https://iweb.dl.sourceforge.net/project/mingw/MinGW/Base/gmp/gmp-5.1.2/gmp-5.1.2-1-mingw32-dll.tar.lzma', '--output', './MinGW/gmp-5.1.2-1-mingw32-dll.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', 'https://managedway.dl.sourceforge.net/project/mingw/MinGW/Base/pthreads-w32/pthreads-w32-2.9.1/pthreads-w32-2.9.1-1-mingw32-dev.tar.lzma', '--output', './MinGW/pthreads-w32-2.9.1-1-mingw32-dev.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://phoenixnap.dl.sourceforge.net/project/mingw/MinGW/Base/pthreads-w32/pthreads-w32-2.9.1/pthreads-w32-2.9.1-1-mingw32-dll.tar.lzma", '--output', './MinGW/pthreads-w32-2.9.1-1-mingw32-dll.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', 'https://pilotfiber.dl.sourceforge.net/project/mingw/MinGW/Base/libiconv/libiconv-1.14-3/libiconv-1.14-3-mingw32-dll.tar.lzma', '--output', './MinGW/libiconv-1.14-3-mingw32-dll.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://gigenet.dl.sourceforge.net/project/mingw/MinGW/Base/zlib/zlib-1.2.8/zlib-1.2.8-1-mingw32-dll.tar.lzma", '--output', './MinGW/zlib-1.2.8-1-mingw32-dll.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://netactuate.dl.sourceforge.net/project/mingw/MinGW/Base/gettext/gettext-0.18.3.1-1/gettext-0.18.3.1-1-mingw32-dll.tar.lzma", '--output', './MinGW/gettext-0.18.3.1-1-mingw32-dll.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://phoenixnap.dl.sourceforge.net/project/mingw/MinGW/Base/gcc/Version4/gcc-4.8.1-4/gcc-core-4.8.1-4-mingw32-bin.tar.lzma", '--output', './MinGW/gcc-core-4.8.1-4-mingw32-bin.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://pilotfiber.dl.sourceforge.net/project/mingw/MinGW/Base/gcc/Version4/gcc-4.8.1-4/gcc-core-4.8.1-4-mingw32-dev.tar.lzma", '--output', './MinGW/gcc-core-4.8.1-4-mingw32-dev.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://netactuate.dl.sourceforge.net/project/mingw/MinGW/Base/gcc/Version4/gcc-4.8.1-4/gcc-core-4.8.1-4-mingw32-dll.tar.lzma", '--output', './MinGW/gcc-core-4.8.1-4-mingw32-dll.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://master.dl.sourceforge.net/project/mingw/MinGW/Base/gcc/Version4/gcc-4.8.1-4/gcc-fortran-4.8.1-4-mingw32-bin.tar.lzma", '--output', './MinGW/gcc-fortran-4.8.1-4-mingw32-bin.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://master.dl.sourceforge.net/project/mingw/MinGW/Base/gcc/Version4/gcc-4.8.1-4/gcc-fortran-4.8.1-4-mingw32-dev.tar.lzma", '--output', './MinGW/gcc-fortran-4.8.1-4-mingw32-dev.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', "https://master.dl.sourceforge.net/project/mingw/MinGW/Base/gcc/Version4/gcc-4.8.1-4/gcc-fortran-4.8.1-4-mingw32-dll.tar.lzma", '--output', './MinGW/gcc-fortran-4.8.1-4-mingw32-dll.tar.lzma'])
            print('\rFetching MinGW... {}%'.format(progress), end='')
            sys.stdout.flush()
            progress += 5
            sp.call(['curl', '-s', 'https://iweb.dl.sourceforge.net/project/mingw/MinGW/Extension/make/make-3.82.90-cvs/make-3.82.90-2-mingw32-cvs-20120902-bin.tar.lzma', '--output', './MinGW/make-3.82.90-2-mingw32-cvs-20120902-bin.tar.lzma'])
            print('\rFetching MinGW...     finished')
        except:
            print("Unexpected error:", sys.exc_info()[0])


        try:
            filenames_xz = ['binutils-2.28-1-mingw32-bin',  'mingwrt-5.0.2-dev', 'mingwrt-5.0.2-dll' ]

            filenames_lzma = ['gcc-core-4.8.1-4-mingw32-bin', 'gcc-core-4.8.1-4-mingw32-dev', 'gcc-core-4.8.1-4-mingw32-dll', 'gcc-fortran-4.8.1-4-mingw32-bin','gcc-fortran-4.8.1-4-mingw32-dev','gcc-fortran-4.8.1-4-mingw32-dll', 'gettext-0.18.3.1-1-mingw32-dll', 'gmp-5.1.2-1-mingw32-dll', 'libiconv-1.14-3-mingw32-dll', 'make-3.82.90-2-mingw32-cvs-20120902-bin', 'mpc-1.0.1-2-mingw32-dll', 'mpfr-3.1.2-2-mingw32-dll', 'pthreads-w32-2.9.1-1-mingw32-dev', 'pthreads-w32-2.9.1-1-mingw32-dll', 'w32api-3.17-2-mingw32-dev', 'zlib-1.2.8-1-mingw32-dll']


            progress = 0
            f = open('mingw_install_output.txt', 'w')
            f2 = open('mingw_install_err.txt','w')
            print('Installing MinGW... {}%'.format(progress), end='')

            for i in range(len(filenames_xz)):
                sp.call(['7z', 'x', filenames_xz[i] + '.tar.xz'], cwd='MinGW', stdout=f, stderr=f2)
                sp.call(['7z', 'x', filenames_xz[i] + '.tar'], cwd='MinGW', stdout=f)
                progress += 5
                print('\rInstalling MinGW... {}%'.format(progress), end='')
                sys.stdout.flush()

            for i in range(len(filenames_lzma)):
                sp.call(['7z', 'x', filenames_lzma[i] + '.tar.lzma'], cwd='MinGW', stdout=f, stderr=f2)
                sp.call(['7z', 'x', filenames_lzma[i] + '.tar'], cwd='MinGW', stdout=f)
                progress += 5
                print('\rInstalling MinGW... {}%'.format(progress), end='')
                sys.stdout.flush()
            print('\rInstalling MinGW... 100%')
            print('Installing MinGW...     finalizing', end='')
            sys.stdout.flush()
            sp.call(['robocopy', 'MinGW', 'C:\\MinGW', '/E'], shell=True, stdout=f)
            path = os.environ.get('path')
            sp.call(['setx.exe','PATH', '{};C:\\MinGW\\bin'.format(path), '/M'], stdout=f)
            print('\rInstalling MinGW...        finished')
            print("Dont forget to restart your shell & run as administrator to access mingw32-make in the future!" )
            f.close()
            f2.close()

            print('Checking for MinGW32-make again... ', end='')
            f = open('check.txt', 'w')
            path  = os.environ.get('PATH')
            sp.call(['mingw32-make'], env=dict(os.environ, PATH='{};C:\\Program Files\\7-Zip;C:\\MinGW\\bin'.format(path)), stdout=f)
            f.close()
            print('    finished')
        except:
            print("Unexpected error:", sys.exc_info()[0])



    #fetch CPT16.5.8 from simon's repository on iri.columbia.edu
    try:
        print('Fetching CPT...', end=' ') #prints Fetching CPT...  without moving to a new line
        sys.stdout.flush() #Forces above line to print - it normally wouldnt until 'finished' prints
        f= open('delouputu.txt', 'w')
        sp.call(['del', '/s', '/q', 'CPT.16.5.8.tar.gz'],stdout=f, shell=True) #deletes CPT file if you already have it
        f.close()
        sp.call(['curl','-s', 'https://iri.columbia.edu/~simon/CPT/CPT.16.5.8.tar.gz', '--output', './CPT.16.5.8.tar.gz']) #curl is a program that grabs data from a provided website url. -s mutes its output. --output specifies the output file.
        print('     Finished') #prints 'Finished'
    except:
        print("Unexpected error:", sys.exc_info()[0]) #if an error occurs, prints the type of error, helpful for debugging


    try:
        print('Decompressing CPT to ./CPT1658...', end='') #prints Decompressing CPT to ./CPT1658...   without moving to a new line
        f = open('cpt_install_output.txt', 'w')
        sys.stdout.flush()  #Forces above line to print - it normally wouldnt until 'finished' prints
        sp.call(['del', '/s', '/q', 'CPT1658'],stdout=f, stderr=f, shell=True) #deletes the directory ./CPT1658 if it exists - just so nothing gets messed up
        sp.call(['rmdir', '/s', '/q', 'CPT1658'],stdout=f, stderr=f, shell=True)
        sp.call(['mkdir', 'CPT1658'], shell=True) # makes directory ./CPT1658 so we can send CPT files there
        f2 = open('./cwd.txt', 'w')
        sp.call(['echo', '%cd%'], stdout=f2, shell=True)
        f2.close()
        f2 = open('./cwd.txt', 'r')
        cwd = f2.readline().strip()
        cwd = cwd.replace('\\', '/')
        f2.close()
        if had_to_get_sevenzip:
            path = os.environ.get('PATH')
            sp.call(['7z','e', 'CPT.16.5.8.tar.gz', '-y'], env=dict(os.environ, PATH='{};C:\\Program Files\\7-Zip'), stdout=f, stderr=f) #Decompresses the .tar.gz CPT file we just downloaded, sends output to direcotry we just made
            sp.call(['7z','x', 'CPT.16.5.8.tar', '-y', '-o{}/CPT1658'.format(cwd)], env=dict(os.environ, PATH='{};C:\\Program Files\\7-Zip'), stdout=f, stderr=f,shell=True) #Decompresses the .tar.gz CPT file we just downloaded, sends output to direcotry we just made
        else:
            path = os.environ.get('PATH')
            sp.call(['7z','e', 'CPT.16.5.8.tar.gz', '-y'], env=dict(os.environ, PATH='{};C:\\Program Files\\7-Zip'), stdout=f, stderr=f) #Decompresses the .tar.gz CPT file we just downloaded, sends output to direcotry we just made
            sp.call(['7z','x', 'CPT.16.5.8.tar','-y', '-o{}/CPT1658/'.format(cwd)], env=dict(os.environ, PATH='{};C:\\Program Files\\7-Zip'), stdout=f, stderr=f,shell=True) #Decompresses the .tar.gz CPT file we just downloaded, sends output to direcotry we just made
        f.close()
        print('     Finished')  #prints 'Finished'
    except:
        print("Unexpected error:", sys.exc_info()[0])


    try:
        comparison_lines = ["gfortran -O2 -frecursive -c -o sggev.o sggev.f","gfortran -O2 -frecursive -c -o sorgtr.o sorgtr.f","gfortran -O2 -frecursive -c -o ssytrs_rook.o ssytrs_rook.f","gfortran -O2 -frecursive -c -o ssyev_2stage.o ssyev_2stage.f","gfortran -O2 -frecursive -c -o dlansy.o dlansy.f","gfortran -O2 -frecursive -c -o dptcon.o dptcon.f","gfortran -O2 -frecursive -c -o dtrti2.o dtrti2.f","gfortran -O2 -frecursive -c -o cgelqf.o cgelqf.f","gfortran -O2 -frecursive -c -o chpevx.o chpevx.f","gfortran -O2 -frecursive -c -o clatrz.o clatrz.f","gfortran -O2 -frecursive -c -o ctrexc.o ctrexc.f","gfortran -O2 -frecursive -c -o zgbbrd.o zgbbrd.f","gfortran -O2 -frecursive -c -o zhetf2_rook.o zhetf2_rook.f","gfortran -O2 -frecursive -c -o zlarfx.o zlarfx.f","gfortran -O2 -frecursive -c -o zsytrf_rk.o zsytrf_rk.f","gfortran -O2 -frecursive -c -o ztprfb.o ztprfb.f","gfortran -O2 -frecursive -c -o dbdsdc.o dbdsdc.f","gfortran -O2 -frecursive -c -o srotmg.o srotmg.f","gfortran -O2 -frecursive -c -o zdscal.o zdscal.f","gfortran -c -O -DDP=1 -DGFORTRAN  -std=f2008 -fall-intrinsics pfv.F95"] #this is for progress checking
        progress = 0
        path = os.environ.get('PATH')
        print('Compiling CPT... {}%'.format(progress), end='') #prints COmpiling CPT
        sys.stdout.flush() #forces above line to print, it wouldn't normally
        f = open('mingw32-make_output.txt', 'w')
        #sp.call(['mingw32-make'], cwd='CPT1658/CPT/16.5.8/lapack/lapack', env=dict(os.environ, PATH='{};C:\\MinGW\\bin'), stdout=f, stderr=f) #remakes lapack - fixes up compile errors
        path = os.environ.get('PATH')
        pipo = sp.Popen(['mingw32-make'], cwd='CPT1658/CPT/16.5.8/lapack/lapack', env=dict(os.environ, PATH='{};C:\\MinGW\\bin'.format(path)), stdout=sp.PIPE, universal_newlines=True) #remakes lapack - fixes up compile errors)
        for line in iter(pipo.stdout.readline, ""): #checks if each line of the output is in the 'Comparison lines' list above- they are 5% increments of progress
            if line.strip() in comparison_lines:
                progress += 5
                print('\rCompiling CPT... {}%'.format(progress), end='') #prints COmpiling CPT
                sys.stdout.flush() #forces above line to print, it wouldn't normally
        pipo.communicate()
        path = os.environ.get('PATH')
        pipo = sp.Popen(['mingw32-make'], cwd='CPT1658/CPT/16.5.8/', env=dict(os.environ, PATH='{};C:\\MinGW\\bin'.format(path)), stdout=sp.PIPE, universal_newlines=True) #compiles CPT using make'
        for line in iter(pipo.stdout.readline, ""):
            if line.strip() in comparison_lines: #if the output line is in the predeterminied group, add 5% to progress (every 105/2100 lines of output this happens)
                progress += 5
                print('\rCompiling CPT... {}%'.format(progress), end='') #prints COmpiling CPT
                sys.stdout.flush() #forces above line to print, it wouldn't normally
        pipo.communicate()
        f.close()

        print('\rCompiling CPT...     Finished') #prints Finished
        f2 = open('./cwd.txt', 'w')
        sp.call(['echo', '%cd%'], stdout=f2, shell=True)
        f2.close()
        f2 = open('./cwd.txt', 'r')
        cwd = f2.readline().strip()
        cwd = cwd.replace('\\', '/')
        f2.close()
        print('Your "cptdir" is {}/CPT1658/CPT/16.5.8/'.format(cwd))
    except sp.CalledProcessError as e:
        print(e.output)
    except:
        print("Unexpected error:", sys.exc_info()[0])



    try:
        print('Checking for conda install... ', end='')
        sys.stdout.flush()  #Forces above line to print - it normally wouldnt until 'finished' prints

        pipo = sp.Popen(['conda'], stdout=sp.PIPE)
        output = pipo.communicate()
        print('     finished')
    #except sp.CalledProcessError as e:
    #    print(e.output)
    except:
        print("Unexpected error:", sys.exc_info()[0])

        try:
            print('Checking for conda install somewhere else... ', end='')
            sys.stdout.flush()
            path = os.environ.get('PATH')
            pipo = sp.Popen(['conda'], env=dict(os.environ, PATH='{};C:/miniconda/bin'.format(path)), stdout=sp.PIPE)
            output = pipo.communicate()
            print('     finished')
        #except sp.CalledProcessError as e:
        #    print(e.output)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print('Installing Miniconda for you... ', end='')
            sys.stdout.flush()
            #f = open(os.path.expanduser('~/miniconda.sh'),'w')
            sp.call(['curl', '-s', 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86.exe', '--output', './Miniconda3-latest-Windows-x86.exe'])#, stdout=f)
            #f.close()
            f = open('./miniconda_install_output.txt', 'w')
            sp.call(['robocopy', '.', 'C:/', 'Miniconda3-latest-Windows-x86.exe'], shell=True,stdout=f)
            sp.call(['Miniconda3-latest-Windows-x86.exe', '/S', '/InstallationType=JustMe', '/D=%UserProfile%\\Miniconda3' ], shell=True, cwd='C:/', stdout=f)
            print('\rInstalling Miniconda for you...      finished')
            print('To access conda, jupyter etc, restart your terminal once this script finishes!')

            sp.call(['setx',  'PATH', '%PATH%;%UserProfile%\\Miniconda3\\Scripts;%UserProfile%\\Miniconda3;%UserProfile%\\Miniconda3\\Library\\bin', '/M'], stdout=f)
            f3 = open('./userprofile.txt', 'w')
            sp.call(['echo', '%UserProfile%'],shell=True, stdout=f3)
            f3.close()
            f3 = open('./userprofile.txt', 'r')
            up = f3.readline().strip()
            up2 = up + '\\Miniconda3\\Library\\bin\\'
            f3.close()
            sp.call(['robocopy', '.', '../../DLLs/', 'libcrypto-1_1.*', 'libssl-1_1.*'], stdout=f, shell=True, cwd=up2)

            print('Checking for conda again... ', end='')
            sys.stdout.flush()
            try:
                f = open('./garbage.txt', 'w')
                sp.call(['conda'], env=dict(os.environ, PATH='{};{}/Miniconda3/Scripts;{}/Miniconda3;{}/Miniconda3/Library/bin'.format(path,up, up, up)),shell=True, stdout=f)
                f.close()
                print('you did it lol')
            except:
                print("Unexpected error:", sys.exc_info()[0])



    try:
        f3 = open('./userprofile.txt', 'w')
        sp.call(['echo', '%UserProfile%'],shell=True, stdout=f3)
        f3.close()
        f3 = open('./userprofile.txt', 'r')
        up = f3.readline().strip()
        up2 = up + '\\Miniconda3\\Library\\bin\\'
        f3.close()
        print('Installing PyCPT Dependencies... ', end='')
        sys.stdout.flush()
        path = os.environ.get('PATH')
        f = open('./conda_output.txt', 'w')
        sp.call(['conda', 'install',  '-y', '-q', 'xarray'],env=dict(os.environ,  PATH='{};{}/Miniconda3/Scripts;{}/Miniconda3'.format(path, up, up)),shell=True, stdout=f, stderr=f)
        print('\rInstalling PyCPT Dependencies... 33%', end='')
        sp.call(['conda', 'install', '-y', '-q', 'jupyter'],env=dict(os.environ, PATH='{};{}/Miniconda3/Scripts;{}/Miniconda3'.format(path, up, up)), shell=True, stdout=f, stderr=f)
        print('\rInstalling PyCPT Dependencies... 66%', end='')
        sp.call(['conda', 'install', '-y', '-q', '-c', 'conda-forge', 'cartopy'], env=dict(os.environ, PATH='{};{}/Miniconda3/Scripts;{}/Miniconda3'.format(path, up, up)),shell=True,stderr=f, stdout=f)
        f.close()
        print('\rInstalling PyCPT Dependencies...      finished')
    except sp.CalledProcessError as e:
        print(e.output)
    except:
        print("Unexpected error:", sys.exc_info()[0])

    try:
        print('Downloading PyCPT... ', end='')
        sp.call(['curl','-s', 'https://raw.githubusercontent.com/kjhall01/PyCPT-Dev/master/scripts/pycpt_functions_seasonal.py', '--output', './pycpt_functions_seasonal.py']) #downloads pycpt_functions_seasonal.py
        sp.call(['curl','-s', 'https://raw.githubusercontent.com/kjhall01/PyCPT-Dev/master/scripts/PyCPT_seav1.4.ipynb', '--output', './PyCPT_seav1.4.ipynb']) #curl is a program that grabs data from a provided website url. -s mutes its output. --output specifies the output file.
        print('      Finished')
    except sp.CalledProcessError as e:
        print(e.output)
    except:
        print("Unexpected error:", sys.exc_info()[0])




























else:
    #fetch CPT16.5.8 from simon's repository on iri.columbia.edu
    try:
        print('Fetching CPT...', end=' ') #prints Fetching CPT...  without moving to a new line
        sys.stdout.flush() #Forces above line to print - it normally wouldnt until 'finished' prints
        sp.call(['rm', '-rf', 'CPT.16.5.8.tar.gz']) #deletes CPT file if you already have it
        sp.call(['curl','-s', 'https://iri.columbia.edu/~simon/CPT/CPT.16.5.8.tar.gz', '--output', './CPT.16.5.8.tar.gz']) #curl is a program that grabs data from a provided website url. -s mutes its output. --output specifies the output file.
        print('     Finished') #prints 'Finished'
    except:
        print("Unexpected error:", sys.exc_info()[0]) #if an error occurs, prints the type of error, helpful for debugging


    try:
        print('Decompressing CPT to ./CPT1658...', end='') #prints Decompressing CPT to ./CPT1658...   without moving to a new line
        sys.stdout.flush()  #Forces above line to print - it normally wouldnt until 'finished' prints
        sp.call(['rm', '-rf', 'CPT1658']) #deletes the directory ./CPT1658 if it exists - just so nothing gets messed up
        sp.call(['mkdir', 'CPT1658']) # makes directory ./CPT1658 so we can send CPT files there
        sp.call(['tar','xf', 'CPT.16.5.8.tar.gz', '--directory', 'CPT1658']) #Decompresses the .tar.gz CPT file we just downloaded, sends output to direcotry we just made
        print('     Finished')  #prints 'Finished'
    except:
        print("Unexpected error:", sys.exc_info()[0])


    try:
        comparison_lines = ["gfortran -O2 -frecursive -c -o sggev.o sggev.f","gfortran -O2 -frecursive -c -o sorgtr.o sorgtr.f","gfortran -O2 -frecursive -c -o ssytrs_rook.o ssytrs_rook.f","gfortran -O2 -frecursive -c -o ssyev_2stage.o ssyev_2stage.f","gfortran -O2 -frecursive -c -o dlansy.o dlansy.f","gfortran -O2 -frecursive -c -o dptcon.o dptcon.f","gfortran -O2 -frecursive -c -o dtrti2.o dtrti2.f","gfortran -O2 -frecursive -c -o cgelqf.o cgelqf.f","gfortran -O2 -frecursive -c -o chpevx.o chpevx.f","gfortran -O2 -frecursive -c -o clatrz.o clatrz.f","gfortran -O2 -frecursive -c -o ctrexc.o ctrexc.f","gfortran -O2 -frecursive -c -o zgbbrd.o zgbbrd.f","gfortran -O2 -frecursive -c -o zhetf2_rook.o zhetf2_rook.f","gfortran -O2 -frecursive -c -o zlarfx.o zlarfx.f","gfortran -O2 -frecursive -c -o zsytrf_rk.o zsytrf_rk.f","gfortran -O2 -frecursive -c -o ztprfb.o ztprfb.f","gfortran -O2 -frecursive -c -o dbdsdc.o dbdsdc.f","gfortran -O2 -frecursive -c -o srotmg.o srotmg.f","gfortran -O2 -frecursive -c -o zdscal.o zdscal.f","gfortran -c -O -DDP=1 -DGFORTRAN  -std=f2008 -fall-intrinsics pfv.F95"] #this is for progress checking
        progress = 0

        print('Compiling CPT... {}%'.format(progress), end='') #prints COmpiling CPT
        sys.stdout.flush() #forces above line to print, it wouldn't normally

        sp.call(['make', '-s', 'clean'], cwd='CPT1658/CPT/16.5.8/lapack/lapack') #remakes lapack - fixes up compile errors
        pipo = sp.Popen(['make'], cwd='CPT1658/CPT/16.5.8/lapack/lapack', stdout=sp.PIPE, universal_newlines=True) #remakes lapack - fixes up compile errors)
        for line in iter(pipo.stdout.readline, ""): #checks if each line of the output is in the 'Comparison lines' list above- they are 5% increments of progress
            if line.strip() in comparison_lines:
                progress += 5
                print('\rCompiling CPT... {}%'.format(progress), end='') #prints COmpiling CPT
                sys.stdout.flush() #forces above line to print, it wouldn't normally
        pipo.communicate()
        pipo = sp.Popen(['make'], cwd='CPT1658/CPT/16.5.8/', stdout=sp.PIPE, universal_newlines=True) #compiles CPT using make'
        for line in iter(pipo.stdout.readline, ""):
            if line.strip() in comparison_lines: #if the output line is in the predeterminied group, add 5% to progress (every 105/2100 lines of output this happens)
                progress += 5
                print('\rCompiling CPT... {}%'.format(progress), end='') #prints COmpiling CPT
                sys.stdout.flush() #forces above line to print, it wouldn't normally
        pipo.communicate()
        print('\rCompiling CPT...     Finished') #prints Finished
        f = open('cptdir.txt', 'w')
        sp.call(['pwd'], stdout=f)
        f.close()
        f = open('cptdir.txt', 'r')
        print('Your "cptdir" is {}CPT1658/CPT/16.5.8/'.format(f.read().strip()))
    except sp.CalledProcessError as e:
        print(e.output)
    except:
        print("Unexpected error:", sys.exc_info()[0])

    try:
        print('Checking for conda install... ', end='')
        sys.stdout.flush()  #Forces above line to print - it normally wouldnt until 'finished' prints

        pipo = sp.Popen(['conda'], stdout=sp.PIPE)
        output = pipo.communicate()
        print('     finished')
    except sp.CalledProcessError as e:
        print(e.output)
    except:
        print("Unexpected error:", sys.exc_info()[0])
        try:
            print('Checking for conda install somewhere else... ', end='')
            sys.stdout.flush()
            os.environ['PATH'] = os.path.expanduser('~/miniconda/bin:{}'.format(os.environ['PATH']))
            pipo = sp.Popen(['conda'], stdout=sp.PIPE)
            output = pipo.communicate()
            print('     finished')
        except sp.CalledProcessError as e:
            print(e.output)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print('Installing Miniconda for you... ')
            sys.stdout.flush()
            f = open(os.path.expanduser('~/miniconda.sh'),'w')
            sp.call(['curl', '-s', 'https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh', '--output', os.path.expanduser('~/miniconda.sh')])#, stdout=f)
            f.close()
            f = open('./miniconda_install_output.txt', 'w')
            sp.call(['bash', os.path.expanduser('~/miniconda.sh'), '-b', '-p', os.path.expanduser('~/miniconda')],stdout=f)
            f.close()
            sys.stdout.write("\033[F")
            print('Installing Miniconda for you...      finished')
            print('Editing your .bash_profile...      finished')
            print('To access conda, jupyter etc, restart your terminal once this script finishes!')
            f = open(os.path.expanduser('~/.bash_profile'), 'a')
            f.write('export PATH=\"{}/miniconda/bin:$PATH\"'.format(os.path.expanduser('~')))
            f.close()
            print('Checking for conda again... ', end='')
            sys.stdout.flush()
            try:
                f = open('./garbage.txt', 'w')
                sp.call(['conda'], stdout=f)
                f.close()
                print('you did it lol')
            except:
                print("Unexpected error:", sys.exc_info()[0])



    try:
        print('Installing PyCPT Dependencies... ', end='')
        sys.stdout.flush()
        f = open('./conda_output.txt', 'w')
        sp.call(['conda', 'install',  '-y', '-q', 'xarray'], stdout=f)
        print('\rInstalling PyCPT Dependencies... 33%', end='')
        sp.call(['conda', 'install', '-y', '-q', 'jupyter'], stdout=f)
        print('\rInstalling PyCPT Dependencies... 66%', end='')
        sp.call(['conda', 'install', '-y', '-q', '-c', 'conda-forge', 'cartopy'], stdout=f)
        f.close()
        print('\rInstalling PyCPT Dependencies...      finished', end='')
    except sp.CalledProcessError as e:
        print(e.output)
    except:
        print("Unexpected error:", sys.exc_info()[0])

    try:
        print('Downloading PyCPT... ', end='')
        sp.call(['curl','-s', 'https://raw.githubusercontent.com/kjhall01/PyCPT-Dev/master/scripts/pycpt_functions_seasonal.py', '--output', './pycpt_functions_seasonal.py']) #downloads pycpt_functions_seasonal.py
        sp.call(['curl','-s', 'https://raw.githubusercontent.com/kjhall01/PyCPT-Dev/master/scripts/PyCPT_seav1.4.ipynb', '--output', './PyCPT_seav1.4.ipynb']) #curl is a program that grabs data from a provided website url. -s mutes its output. --output specifies the output file.
        print('      Finished')
    except sp.CalledProcessError as e:
        print(e.output)
    except:
        print("Unexpected error:", sys.exc_info()[0])
