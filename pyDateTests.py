
import os;
import pyDate;


def testPyDate(tyear,tmonth,tday,tdoy,tgpsWeek,tgpsWeekDay):
    
    print

    # convert year month day to day of year
    doy = pyDate.date2doy(tyear,tmonth,tday)[0]
    
    # test + write
    if doy != tdoy:
        os.sys.stderr.write("ymd to doy failed!\n");
    else:
        os.sys.stdout.write('ymd to day passed!\n');
    
    # print the result
    print tyear,tmonth,tday,"-->", tyear,doy
    print
    
    # convert year and day of year to year month day
    month,day = pyDate.doy2date(tyear, doy);
    
    # test
    if month != tmonth or day != tday:
        os.sys.stderr.write("year doy to year month day failed!\n");
    else:
        os.sys.stdout.write("year doy to year month day passed!\n")
    
    # print the results
    print tyear,doy ,"-->", tyear,month,day
    print 
    
    gpsWeek,gpsWeekDay = pyDate.date2gpsDate(tyear,tmonth,tday);
    
    if gpsWeek != tgpsWeek or gpsWeekDay != tgpsWeekDay:
        os.sys.stderr.write("year, month,day to gps date failed!\n");
    else:
        os.sys.stdout.write("year, month, day to gps date passed!\n");
        
    print tyear,tmonth,tday,"-->", gpsWeek, gpsWeekDay
    print
    
    mjd = pyDate.gpsDate2mjd(gpsWeek,gpsWeekDay);
    year,month,day = pyDate.mjd2date(mjd);
    
    if year != tyear or day != tday or month != tmonth:
        os.sys.stderr.write("gps date to date failed!\n");
    else:
        os.sys.stdout.write("gps date to date passed!\n");
        
    print gpsWeek, gpsWeekDay, "-->", year, month, day;
    print
    
    # ok, now test object initialization
    pydate = pyDate.Date(year=year,doy=doy);
    if pydate.year != tyear or pydate.doy != tdoy \
        or pydate.month != tmonth or pydate.day !=tday \
            or pydate.gpsWeek !=tgpsWeek or pydate.gpsWeekDay != tgpsWeekDay:
        os.sys.stderr.write('pydate year,doy initialization failed\n')
    else:
        os.sys.stdout.write('pyDate.Date() year,doy init passed\n');
        
    pydate = pyDate.Date(year=year,day=day,month=month);
    if pydate.year != tyear or pydate.doy != tdoy \
        or pydate.month != tmonth or pydate.day !=tday \
            or pydate.gpsWeek !=tgpsWeek or pydate.gpsWeekDay != tgpsWeekDay:
        os.sys.stderr.write('pydate year,month,day initialization failed\n')
    else:
        os.sys.stdout.write('pyDate.Date() year,month,day init passed\n');
        
    pydate = pyDate.Date(gpsWeek=gpsWeek,gpsWeekDay=gpsWeekDay);
    if pydate.year != tyear or pydate.doy != tdoy \
        or pydate.month != tmonth or pydate.day !=tday \
            or pydate.gpsWeek !=tgpsWeek or pydate.gpsWeekDay != tgpsWeekDay:
        os.sys.stderr.write('pydate gpsWeek,gpsWeekDay initialization failed\n')
    else:
        os.sys.stdout.write('pyDate.Date() gpsWeek,gpsWeekDay init passed\n');


    # test operator overloads and rich comparison operators
    print
    if (pyDate.Date(year=tyear,doy=tdoy) < pyDate.Date(year=tyear,doy=tdoy)-1) == True:
        os.sys.stderr.write('pyDate.Date lt operator failed!\n');
    else:
        os.sys.stdout.write('pyDate.Date lt operator passed!\n');

    if (pyDate.Date(year=tyear,doy=tdoy) > pyDate.Date(year=tyear,doy=tdoy)-1) == False:
        os.sys.stderr.write('pyDate.Date gt operator failed!\n');
    else:
        os.sys.stdout.write('pyDate.Date gt operator passed!\n');
        
    if (pyDate.Date(year=tyear,doy=tdoy) <= pyDate.Date(year=tyear,doy=tdoy)-1) == True:
        os.sys.stderr.write('pyDate.Date le operator failed!\n');
    else:
        os.sys.stdout.write('pyDate.Date le operator passed!\n');
        
    if (pyDate.Date(year=tyear,doy=tdoy) >= pyDate.Date(year=tyear,doy=tdoy)-1) == False:
        os.sys.stderr.write('pyDate.Date ge operator failed!\n');
    else:
        os.sys.stdout.write('pyDate.Date ge operator passed!\n');
        
    if (pyDate.Date(year=tyear,doy=tdoy) == pyDate.Date(year=tyear,doy=tdoy)) == False:
        os.sys.stderr.write('pyDate.Date eq operator failed!\n');
    else:
        os.sys.stdout.write('pyDate.Date eq operator passed!\n');
        
    if (pyDate.Date(year=tyear,doy=tdoy) != pyDate.Date(year=tyear,doy=tdoy)) == True:
        os.sys.stderr.write('pyDate.Date ne operator failed!\n');
    else:
        os.sys.stdout.write('pyDate.Date ne operator passed!\n');
        
        
testPyDate(2006,2,26,57,1364,0);
testPyDate(2004,8,31,244,1286,2);
testPyDate(1994,9,12,255,766,1);

testPyDate(2010,12,31,365,1616,5);
testPyDate(2011,1,1,1,1616,6);

testPyDate(1996,12,31,366,886,2);
testPyDate(1997,1,1,1,886,3);

print pyDate.Date(year=2005,doy=1).fyear


