#include <cmath>
#include <fstream> // this includes iostream.h
#include <iomanip>

#include "consts.h"
#include "datetime.h"
#include "intrpsp3.h"
#include "physcon.h"
#include "rinex.h"

using std::string;
using namespace NGSdatetime;
using namespace NGSrinex;

void llh2xyz(double a, double finverse, double lat, double lon, double h,
             double *x, double *y, double *z);

void xyz2llh(double a, double finverse, double x, double y, double z,
             double *lat, double *lon, double *h);

double dot(int n, double a[], double b[]);

int bccalc(double tsubr, int isv, double recf[3], double vecf[3]);
int bcorb(double tc, int isv, double recf[4], double vecf[3]);
int bcread(RinexNavFile &navFile, ofstream &out);

string padZeros(int a);

const int MAXSVS = 36;
const int MAXUNK = 28;
const int MAXEPOCH = 50; // PRN blocks per file, was 160

double xco[MAXSVS][MAXUNK][MAXEPOCH];
int nxco[MAXSVS];
int jxco;
long cuml_bad_blocks[MAXSVS];
double bcpos[4], bcvel[3];
double t_r;
int BRprns[MAXSVS];
#pragma argsused
int main() {
  int i, j, ierr, prnNum;
  double jpi = 4.0 * atan(1.0);
  SP3File mysp3;
  double tt;
  double lat, lon, h;
  double pvVec[10];
  double svpos[10];
  double n[3], e[3], u[3], rho[3], len, zenithAng, tmpe, tmpn, azim, elvAng;
  double xsta[3];
  string orbfile, outfile, stemp, rinfile;
  YMDHMS ymdhms;
  double tmp1, tmp2, sz;
  double xmap, ymap;
  double r[20];
  YMDHMS prevYmdhms;
  double year1, mon1, day1, hr1, min1, sec1;
  double year2, mon2, day2, hr2, min2, sec2;
  DateTime startTime;
  DateTime endTime;
  DateTime currEpoch;
  bool broadcastExist = false;
  double cutoffAngle;
  double arrowsx[40][2];
  double arrowsy[40][2];
  string tString;
  long index;
  GPSTime gpsTime;
  double htt;
  double oneMinute = 5 / 60. / 24.;
  int icnt = 0;
  double dist3d;
  bool LLHsupplied;
  double xtmp, ytmp, ztmp;

  for (i = 0; i < 40; i++) {
    arrowsx[i][0] = -9999.0;
    arrowsx[i][1] = -9999.0;
    arrowsy[i][0] = -9999.0;
    arrowsy[i][1] = -9999.0;
  }

  for (i = 1; i < MAXSVS; i++) {
    BRprns[i] = i;
  }

  // ------------------------------------------------------
  // open log file
  // ------------------------------------------------------

  ofstream out("skyplot.log");
  if (!out) {
    cerr << "Error opening skyplot.out ! " << endl;
    return -1;
  }

  out.setf(ios::fixed, ios::floatfield);

  ifstream inp("skyplot.inp");
  if (!inp) {
    cerr << "Error opening skyplot.inp ! " << endl;
    out << "Error opening skyplot.inp ! " << endl << endl;
    out << "\n\nABNORMAL Termination " << endl;
    out.close();
    return -1;
  }

  inp >> year1 >> mon1 >> day1 >> hr1 >> min1 >> sec1;
  inp >> year2 >> mon2 >> day2 >> hr2 >> min2 >> sec2;
  inp >> orbfile;
  inp >> stemp >> xsta[0] >> xsta[1] >> xsta[2];
  inp >> cutoffAngle;
  inp.close();

  out << "Input file: " << endl;
  out.precision(0);
  out << year1 << " " << mon1 << " " << day1 << " " << hr1 << " " << min1 << " "
      << sec1 << endl;
  out << year2 << " " << mon2 << " " << day2 << " " << hr2 << " " << min2 << " "
      << sec2 << endl;
  out << orbfile << endl;
  out.precision(3);
  out << stemp << " " << xsta[0] << " " << xsta[1] << " " << xsta[2] << " "
      << endl;
  out << cutoffAngle << endl;

  dist3d = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1] + xsta[2] * xsta[2]);

  if (dist3d >= 6300000.0) {
    //  xyz given as input
    xyz2llh(6378137.0, 298.257223563, xsta[0], xsta[1], xsta[2], &lat, &lon,
            &htt);
    out << "\n\nuser supplied XYZ " << endl;
    LLHsupplied = false;
  } else {
    //  llh given as input
    lat = xsta[0] * jpi / 180.0;
    lon = xsta[1] * jpi / 180.0;
    htt = xsta[2];
    out << "\n\nuser supplied LLH " << endl;
    LLHsupplied = true;
  }

  out.precision(7);

  // out << "LLH:  " << lat << " " << lon << " " << htt << endl;

  // ------------------------------------------------------
  // create title for GMT
  // ------------------------------------------------------

  ofstream outtt("title.txt");
  if (!outtt) {
    cerr << "Error opening title.txt ! " << endl;
    return -1;
  }
  outtt.setf(ios::fixed, ios::floatfield);

  outtt << setw(2) << setprecision(0);
  outtt << "0 2.4  12 0 0 CM " << stemp << endl;
  outtt << "0 2.2  7 0 0 CM Lat: " << setw(10) << setprecision(4)
        << lat * 180.0 / jpi << "\312    Lon:  " << setw(10) << setprecision(4)
        << lon * 180.0 / jpi << "\312    Ell Ht:  " << setw(10)
        << setprecision(1) << htt << " (m) " << endl;
  outtt << "0 2.0  7 0 0 CM GPS Time:   Start " << setw(2) << setprecision(0)
        << year1 << "/" << padZeros(mon1) << "/" << padZeros(day1) << "  "
        << padZeros(hr1) << ":" << padZeros(min1) << ":" << padZeros(sec1)
        << "      Stop " << year2 << "/" << padZeros(mon2) << "/"
        << padZeros(day2) << "  " << padZeros(hr2) << ":" << padZeros(min2)
        << ":" << padZeros(sec2) << endl;

  outtt.close();

  // ------------------------------------------------------
  // open the DOS batch file
  // ------------------------------------------------------

  ofstream outbat("skyplot.sh");
  if (!outbat) {
    cerr << "Error opening skyplot.sh ! " << endl;
    return -1;
  }

  outbat << "#!/bin/bash" << endl;
  outbat << "gmt psxy elevRings.dat  -R-1.6/1.6/-1.6/1.6 -JX7.0  -W1.0p,0/0/0 "
            "-G230 -V -K -P -X0.75 -Y1.0 > skyplot.ps "
         << endl;
  outbat
      << "gmt psxy cutoffRing.dat -R -JX  -W0.2p,255/0/0  -G255 -V -O -K -P >> "
         "skyplot.ps "
      << endl;
  outbat << "gmt psxy elevRings.dat  -R -JX  -W1.0p,0/0/0  -V -O -K -P >> "
            "skyplot.ps "
         << endl;
  outbat
      << "gmt psxy elevRings.dat  -R -JX  -W0.5p,255/255/255  -V -O -K -P >> "
         "skyplot.ps "
      << endl;
  outbat << "gmt pstext title.txt   -R -JX  -N -V  -O -K -P >> skyplot.ps "
         << endl;

  // ------------------------------------------------------
  // broadcast or precise?
  // ------------------------------------------------------

  ifstream inp2(orbfile.c_str());
  if (!inp2) {
    cerr << "Warning: " << orbfile << " does not exist " << endl;
    inp2.close();
    exit(0);
  } else {
    getline(inp2, tString);
    index = tString.find("NAV");
    if (0 < index && index < 80) {
      out << "\n\nUsing broadcast file " << orbfile << endl;
      broadcastExist = true;

      RinexNavFile mynav;
      try {
        mynav.setPathFilenameMode(orbfile, ios_base::in);
      } catch (RinexFileException &openExcep) {
        cout << "Error opening file: " << orbfile << endl;
        cout << "Rinex File Exception is: " << endl
             << openExcep.getMessage() << endl;
        cout << endl
             << "Error Messages for " << orbfile << " :" << endl
             << mynav.getErrorMessages() << endl
             << endl;
        cout << endl
             << "Format Warnings for " << orbfile << " :" << endl
             << mynav.getWarningMessages() << endl
             << endl
             << endl;
        cout << "Terminating program POINT due to an error." << endl;
      }

      // cout << endl << "NV::Header for RINEX NAV file: " << orbfile << endl;
      try {
        bcread(mynav, out);
      } catch (RequiredRecordMissingException &headerExcep) {
        cout << "RequiredRecordMissingException is: " << endl
             << headerExcep.getMessage() << endl;
        cout << endl
             << "Error Messages for " << orbfile << " :" << endl
             << mynav.getErrorMessages() << endl
             << endl;
        cout << endl
             << "Format Warnings for " << orbfile << " :" << endl
             << mynav.getWarningMessages() << endl
             << endl
             << endl;
        cout << "Terminating program POINT due to an error." << endl;
      }

    } else {
      cout << "\n\nUsing sp3 file " << orbfile << endl;
      mysp3.setPathFilenameMode(orbfile.c_str(), ios_base::in);
      mysp3.readHeader();
    }
  }
  inp2.close();

  startTime = DateTime(year1, mon1, day1, hr1, min1, sec1);
  endTime = DateTime(year2, mon2, day2, hr2, min2, sec2);

  if (!LLHsupplied) {
    xyz2llh(AE84, FLATINV84, xsta[0], xsta[1], xsta[2], &lat, &lon, &h);
  } else {
    llh2xyz(AE84, FLATINV84, lat, lon, htt, &xtmp, &ytmp, &ztmp);

    xsta[0] = xtmp;
    xsta[1] = ytmp;
    xsta[2] = ztmp;
  }

  n[0] = -sin(lat) * cos(lon);
  n[1] = -sin(lat) * sin(lon);
  n[2] = cos(lat);

  e[0] = -sin(lon);
  e[1] = cos(lon);
  e[2] = 0.0;

  u[0] = cos(lat) * cos(lon);
  u[1] = cos(lat) * sin(lon);
  u[2] = sin(lat);

  cout.setf(ios::fixed, ios::floatfield);

  // ------------------------------------------------------
  // create input files for GMT
  // ------------------------------------------------------

  system(" rm -f *.sat.xy arrows.xy elevRings.dat ring.txt nesw.txt skyplot.ps "
         ".gmt* 2> del.me");
  ofstream satcoords;

  ofstream hourstamps("hr.txt");
  if (!hourstamps) {
    cout << "ERROR: cant open hr.txt " << endl;
    exit(0);
  }

  ofstream hourdots("hr.xy");
  if (!hourdots) {
    cout << "ERROR: cant open hr.xy " << endl;
    exit(0);
  }

  ofstream azelOut("skyplot.AzEl");
  if (!azelOut) {
    cout << "ERROR: cant open skyplot.AzEl " << endl;
    exit(0);
  }

  azelOut.setf(ios::fixed, ios::floatfield);

  // loop over time

  for (currEpoch = startTime; currEpoch <= endTime;
       currEpoch = (currEpoch + oneMinute)) {
    icnt++;

    // integral hour required for plotting

    ymdhms = currEpoch.GetYMDHMS();
    if (ymdhms.hour != prevYmdhms.hour && icnt > 1) {
      currEpoch = DateTime(ymdhms.year, ymdhms.month, ymdhms.day, ymdhms.hour,
                           0, 0.00000001);
      ymdhms = currEpoch.GetYMDHMS();
    }

    gpsTime = currEpoch.GetGPSTime();
    t_r = gpsTime.secsOfWeek; // this is the given receipt time (dt)

    // loop over all satellites

    for (i = 1; i < MAXSVS; i++) {
      prnNum = i;

      if (broadcastExist) {

        ierr = bcorb(t_r, i, bcpos, bcvel); // BR

        if (fabs(bcpos[0]) < 1.0e-8) {
          ierr = 1;
        } else {
          svpos[0] = bcpos[0];
          svpos[1] = bcpos[1];
          svpos[2] = bcpos[2];
          // cout <<  i << "    " << svpos[0] << "  " << svpos[1] << " "
          //  << svpos[2] << endl;
        }
      } else {
        ierr = (int)mysp3.getSVPosVel(currEpoch, (unsigned short)prnNum, pvVec);

        svpos[0] = pvVec[0] * 1000.0;
        svpos[1] = pvVec[1] * 1000.0;
        svpos[2] = pvVec[2] * 1000.0;
      }
      if (ierr == 0) {
        // pg 148 hoffman-wellenhoff

        for (j = 0; j < 3; j++) {
          rho[j] = svpos[j] - xsta[j];
        }
        len = 0.0;
        for (j = 0; j < 3; j++) {
          len += (rho[j] * rho[j]);
        }
        len = sqrt(len);
        for (j = 0; j < 3; j++) {
          rho[j] /= len;
        }
        zenithAng = acos(dot(3, rho, u));

        tmp1 = dot(3, rho, e);
        sz = sin(zenithAng);
        tmpe = tmp1 / sz;
        tmp2 = dot(3, rho, n);
        tmpn = tmp2 / sz;

        azim = atan2(tmpe, tmpn);
        if (azim < 0.0)
          azim += (2.0 * jpi);

        elvAng = jpi / 2.0 - zenithAng;

        if (elvAng * 180.0 / jpi > cutoffAngle) {

          // compute map coords for satellite

          xmap = (jpi / 2.0 - elvAng) * sin(azim);
          ymap = (jpi / 2.0 - elvAng) * cos(azim);

          tt = (ymdhms.hour + ymdhms.min / 60.0 + ymdhms.sec / 3600.0);

          azelOut << setw(10) << setprecision(6) << tt << setw(4)
                  << setprecision(0) << prnNum << setw(9) << setprecision(3)
                  << (azim * 180.0 / jpi) << setw(9) << setprecision(3)
                  << (elvAng * 180.0 / jpi) <<
              // setw(18) << setprecision(13) << xmap <<
              // setw(18) << setprecision(13) << ymap <<
              endl;

          ostringstream fname;
          fname << prnNum << ".sat.xy";

          satcoords.open(fname.str().c_str(), ios::app);
          satcoords << xmap << "    " << ymap << endl;

          satcoords.close();

          if (tt - (double)((int)tt) < 1.0e-6 && ymdhms.min != 59 &&
              ymdhms.sec != 60.000) {

            /* cout << "hourstamp " << tt << " " << currEpoch <<
               "  " << ymdhms.hour <<   endl; */

            hourstamps << (xmap + 0.02) << " " << ymap << " 6 0 0 ML "
                       << (int)ymdhms.hour << "h" << endl;

            hourdots << xmap << " " << ymap << endl;
          }

          // save arrow info
          arrowsx[prnNum][0] = arrowsx[prnNum][1];
          arrowsx[prnNum][1] = xmap;

          arrowsy[prnNum][0] = arrowsy[prnNum][1];
          arrowsy[prnNum][1] = ymap;

        } // end if
      }
    } // end loop over satellites
    prevYmdhms = ymdhms;
  } // end loop over time

  hourstamps.close();
  hourdots.close();
  azelOut.close();

  out << "\n\nFinished writing azimuth/elevation " << endl;

  // ------------------------------------------------------
  // create arrow input files for GMT
  // ------------------------------------------------------

  // arrow data

  ofstream arrowsf("arrows.xy");
  if (!arrowsf) {
    cout << "ERROR: cant open arrows.xy " << endl;
    exit(0);
  }

  double rad = 0.1;
  double dx, dy, az;
  int acount = 0;
  for (i = 0; i < 40; i++) {
    if (arrowsx[i][0] != -9999.0 && arrowsx[i][1] != -9999.0 &&
        arrowsy[i][0] != -9999.0 && arrowsy[i][1] != -9999.0) {
      acount++;
      dx = arrowsx[i][1] - arrowsx[i][0];
      dy = arrowsy[i][1] - arrowsy[i][0];
      az = atan2(dx, dy);
      dx = rad * sin(az);
      dy = rad * cos(az);

      arrowsf << arrowsx[i][1] << "  " << arrowsy[i][1] << " " << dx << " "
              << dy << " 0 0 0 " << i << endl;

      outbat << "gmt psxy " << i
             << ".sat.xy  -R -JX -W0.75p,0/0/0  -V -P -O -K >> skyplot.ps"
             << endl;
    }
  }

  arrowsf.close();

  // ------------------------------------------------------
  // create elevation ring nput files for GMT
  // ------------------------------------------------------

  ofstream outf("elevRings.dat");
  if (!outf) {
    cerr << "Error opening elevRings.dat ! " << endl;
    return -1;
  }

  ofstream outCut("cutoffRing.dat");
  if (!outCut) {
    cerr << "Error opening cutOffRing.dat ! " << endl;
    return -1;
  }

  outf.setf(ios::fixed, ios::floatfield);
  outCut.setf(ios::fixed, ios::floatfield);

  r[0] = jpi / 2.0 - 0. * jpi / 180.0;
  r[1] = jpi / 2.0 - 30. * jpi / 180.0;
  r[2] = jpi / 2.0 - 60. * jpi / 180.0;
  r[3] = jpi / 2.0 - cutoffAngle * jpi / 180.0;

  // nominal rings
  for (j = 0; j < 3; j++) {
    outf << ">" << endl;
    for (i = 0; i <= 360; i += 1) {
      xmap = r[j] * cos(i * jpi / 180.0);
      ymap = r[j] * sin(i * jpi / 180.0);
      outf << setw(18) << setprecision(8) << xmap << setw(18) << setprecision(8)
           << ymap << endl;
    }
  }

  // cutoff rings
  j = 3;
  outf << ">" << endl;
  for (i = 0; i <= 360; i += 1) {
    xmap = r[j] * cos(i * jpi / 180.0);
    ymap = r[j] * sin(i * jpi / 180.0);
    outCut << setw(18) << setprecision(8) << xmap << setw(18) << setprecision(8)
           << ymap << endl;
  }

  outf.close();
  outCut.close();

  ofstream outcross("cross.txt");
  if (!outcross) {
    cerr << "Error opening cross.txt ! " << endl;
    return -1;
  }

  // ------------------------------------------------------
  // create input files for GMT
  // ------------------------------------------------------

  // nesw cross

  outcross << ">" << endl;
  outcross << setw(18) << setprecision(8) << (120 * jpi / 180.0) << setw(18)
           << setprecision(8) << 0 << endl;
  outcross << setw(18) << setprecision(8) << (-120 * jpi / 180.0) << setw(18)
           << setprecision(8) << 0 << endl;

  outcross << ">" << endl;
  outcross << setw(18) << setprecision(8) << 0 << setw(18) << setprecision(8)
           << (120 * jpi / 180.0) << endl;
  outcross << setw(18) << setprecision(8) << 0 << setw(18) << setprecision(8)
           << (-120 * jpi / 180.0) << endl;

  // nw/sw
  double rtemp = 92.0;
  double qtemp = 88.0;
  outcross << ">" << endl;
  outcross << setw(18) << setprecision(8)
           << sin(jpi / 4.0) * (rtemp * jpi / 180.0) << setw(18)
           << setprecision(8) << cos(jpi / 4.0) * (rtemp * jpi / 180.0) << endl;
  outcross << setw(18) << setprecision(8)
           << sin(jpi / 4.0) * (qtemp * jpi / 180.0) << setw(18)
           << setprecision(8) << cos(jpi / 4.0) * (qtemp * jpi / 180.0) << endl;

  outcross << ">" << endl;
  outcross << setw(18) << setprecision(8)
           << -sin(jpi / 4.0) * (rtemp * jpi / 180.0) << setw(18)
           << setprecision(8) << -cos(jpi / 4.0) * (rtemp * jpi / 180.0)
           << endl;
  outcross << setw(18) << setprecision(8)
           << -sin(jpi / 4.0) * (qtemp * jpi / 180.0) << setw(18)
           << setprecision(8) << -cos(jpi / 4.0) * (qtemp * jpi / 180.0)
           << endl;

  outcross << ">" << endl;
  outcross << setw(18) << setprecision(8)
           << -sin(jpi / 4.0) * (rtemp * jpi / 180.0) << setw(18)
           << setprecision(8) << cos(jpi / 4.0) * (rtemp * jpi / 180.0) << endl;
  outcross << setw(18) << setprecision(8)
           << -sin(jpi / 4.0) * (qtemp * jpi / 180.0) << setw(18)
           << setprecision(8) << cos(jpi / 4.0) * (qtemp * jpi / 180.0) << endl;

  outcross << ">" << endl;
  outcross << setw(18) << setprecision(8)
           << sin(jpi / 4.0) * (rtemp * jpi / 180.0) << setw(18)
           << setprecision(8) << -cos(jpi / 4.0) * (rtemp * jpi / 180.0)
           << endl;
  outcross << setw(18) << setprecision(8)
           << sin(jpi / 4.0) * (qtemp * jpi / 180.0) << setw(18)
           << setprecision(8) << -cos(jpi / 4.0) * (qtemp * jpi / 180.0)
           << endl;

  outcross.close();

  // ------------------------------------------------------
  // create input files for GMT
  // ------------------------------------------------------

  // ring labels

  ofstream outring("ring.txt");
  if (!outring) {
    cerr << "Error opening ring.txt ! " << endl;
    return -1;
  }

  double x, y;
  y = 0;
  x = jpi / 2.0 - 0.0;
  //   \312 make the degree symbol
  outring << x << "  " << y << " 8 0 0 CM  0\312" << endl;
  ;

  x = jpi / 2.0 - 30.0 * jpi / 180.0;
  outring << x << "  " << y << " 8 0 0 CM 30\312" << endl;
  ;

  x = jpi / 2.0 - 60.0 * jpi / 180.0;
  outring << x << " " << y << " 8 0 0 CM 60\312" << endl;
  ;

  x = jpi / 2.0 - 90.0 * jpi / 180.0;
  outring << x << " " << y << " 8 0 0 CM 90\312" << endl;
  ;

  x = jpi / 2.0 - cutoffAngle * jpi / 180.0;
  outring << -x << " " << y << " 8 0 0 CM " << cutoffAngle << "\312" << endl;
  ;

  outring.close();

  // NESW labels

  ofstream nesw("nesw.txt");
  if (!nesw) {
    cerr << "Error opening nesw.txt ! " << endl;
    return -1;
  }

  y = 15 * 3.1415 / 180;

  nesw << 0 << " " << 100 * 3.1415 / 180 << "  10 0 0 CM N" << endl;
  nesw << 0 << " " << -100 * 3.1415 / 180 << "  10 0 0 CM S" << endl;
  nesw << -100 * 3.1415 / 180 << " " << 0 << "  10 0 0 CM W" << endl;
  nesw << 100 * 3.1415 / 180 << " " << 0 << "  10 0 0 CM E" << endl;

  rtemp += 4.0;

  x = sin(jpi / 4.0) * (rtemp * jpi / 180.0);
  y = cos(jpi / 4.0) * (rtemp * jpi / 180.0);
  nesw << x << " " << y << "  10 0 0 CM 45\312" << endl;

  x = sin(jpi / 4.0) * (rtemp * jpi / 180.0);
  y = -cos(jpi / 4.0) * (rtemp * jpi / 180.0);
  nesw << x << " " << y << "  10 0 0 CM 135\312" << endl;

  x = -sin(jpi / 4.0) * (rtemp * jpi / 180.0);
  y = -cos(jpi / 4.0) * (rtemp * jpi / 180.0);
  nesw << x << " " << y << "  10 0 0 CM 225\312" << endl;

  x = -sin(jpi / 4.0) * (rtemp * jpi / 180.0);
  y = cos(jpi / 4.0) * (rtemp * jpi / 180.0);
  nesw << x << " " << y << "  10 0 0 CM 315\312" << endl;

  nesw.close();

  outbat
      << "gmt psxy hr.xy  -R -JX -V  -Sc0.03 -G0/0/0  -O -K -P >> skyplot.ps "
      << endl;
  outbat << "gmt pstext hr.txt  -R -JX -V -G0/0/255  -O -K -P >> skyplot.ps "
         << endl;
  outbat << "gmt psxy cross.txt  -R -JX   -V -O -K -P >> skyplot.ps " << endl;
  outbat << "gmt pstext nesw.txt  -R  -JX -O -K  -N   >> skyplot.ps " << endl;
  outbat << "gmt psvelo arrows.xy  -R -JX  -L  -W1.0p,255/0/0 -Se1/0.95/12 "
            "-A0.0020/0.035/0.025 "
         << " -N  -H0 -O -K -P -V >>  skyplot.ps " << endl;

  outbat << "gmt pstext ring.txt  -R  -JX -O   -N -W1p  >> skyplot.ps " << endl;
  outbat << "gmt psconvert skyplot.ps -A -Tg" << endl;

  outbat << "echo ------------------------------ " << endl;
  outbat << "echo ------------------------------ " << endl;
  outbat << "echo \"View or print skyplot.ps (or skyplot.png)\"" << endl;

  outbat.close();

  out << "\n\nNormal Termination " << endl;

  out.close();
  system(" rm -f del.me  "); // cleanup
  // system("skyplot.bat 1> gmt.messages ");
  // system("del .gmt* ring.txt arrows.xy *.sat.xy title.txt hr.txt ");

  cerr << "\n\nNormal Termination.   "
       << "Issue the command 'skyplot.bat' to generate 'skyplot.ps' " << endl;

  return 0;
}

// --------------------------------------------------------------------------

double dot(int n, double a[], double b[]) {
  int i;
  double len;
  len = 0.0;
  for (i = 0; i < n; i++) {
    len += (a[i] * b[i]);
  }
  return len;
}

// --------------------------------------------------------------------------

/* ---------------------------------------------- */
/* Anderson & Mikhail, Surveying Theory & Practice */
/* pg 1124 */

void xyz2llh(double a, double finverse, double x, double y, double z,
             double *lat, double *lon, double *h) {
  int iter;
  double e2, delta, tmp, No, lat1, lat_prev, htmp, f;
  // double pi = 4.0*atan(1.0);
  f = 1.0 / finverse;
  e2 = (2.0 - f) * f;

  /*
    printf("x,y,z %lf %lf %lf\n",x,y,z);
    printf("a,e2 %lf %20.10le\n",a,e2);
  */
  iter = 0;
  delta = 1.0;

  if (sqrt(x * x + y * y) < 1.0e-10) {
    cout << "ERROR: attempting square root of zero in xyz2llh " << endl;
    exit(0);
  }

  lat_prev = atan(z / sqrt(x * x + y * y));
  while (fabs(delta) > 1.0e-12) {
    iter++;
    if (iter > 20) {
      printf("no convergence in xyz2llh\n");
      exit(0);
    }
    No = a / sqrt(1.0 - e2 * sin(lat_prev) * sin(lat_prev));
    htmp = sqrt(x * x + y * y) / cos(lat_prev) - No;
    tmp = 1.0 - e2 * (No / (No + htmp));
    lat1 = atan(z / sqrt(x * x + y * y) / tmp);
    delta = lat_prev - lat1;
    /*
      printf("lat_prev, lat1 %lf %lf %lf\n",
      lat_prev*180.0/pi,lat1*180.0/pi,htmp);
    */
    lat_prev = lat1;
  }
  *lat = lat1;
  *lon = atan2(y, x);
  *h = htmp;
}

//---------------------------------------------------------------------------

/*
c
c********1*********2*********3*********4*********5*********6*********7**
c
c name:        bcorb
c version:     9010.20
c written by:  Dr. G. L. Mader
c purpose:     This subroutine picks the bc elements nearest or just
c              before and after the current time and computes the SV
c              position from each set. A weighted average SV position
c              is then computed.
c
c input parameters
c ------------------------------------
c isv           satellite PRN number
c tc            current time - GPS [seconds of week]
c
c output parameters
c -------------------------------------
c ierr               error flag
c                      not equal to zero = error
c ivel               velocity computation flag
c                      = 1 = compute velocity
c recf(3)            position (meters)
c vecf(3)            velocity (meters/s)
c
c local variables and constants:
c -------------------------------------
c n                 loop counter
c n1                index of element set just before time tc
c n2                index of element set just after time tc
c recf1(3)          xyz coordinates from first set of orbit elements
c recf2(3)          xyz coordinates from second set of orbit elements
c tcx               scratch variable holding tc
c tdif              difference between the epochs of the two
c                   orbit element sets used in computing a linear
c                   weighting factor for averaging the positions
c                   derived from these element sets.
c tk                difference between the requested epoch and the
c                   epoch of the orbit elements used in selecting
c                   the nearest element sets [seconds of week].
c vecf1(3)          xyz velocity from first set of orbit elements
c vecf2(3)          xyz velocity from second set of orbit elements
c w1                weighting factor for element set 1
c w2                weighting factor for element set 2
c xcot              epoch of the orbit element set
c
c global variables and constants:
c  -------------------------------------
c jxco              index of orbit element set selected for computation
c xco()             all sets of broadcast orbit elements
c nxco()            number of sets for each satellite
c
c
c calls:
c -------------------------------------
c bccalc            compute satellite position and velocity from
c                   the broadcast orbit parameters
c
c include files used:
c -------------------------------------
c bcorbit.cm        broadcast orbit and clock elements
c
c common blocks:
c -------------------------------------
c /bcorbit/         broadcast orbit and clock elements
c
c references:
c -------------------------------------
c OMNI Users Guide, April 1989, G.L. Mader, National
c   Geodetic Survey, Rockville, Maryland.
c
c comments:
c -------------------------------------
c
c********1*********2*********3*********4*********5*********6*********7**
c::modification history:
c::8905.03, GLM, Creation of program merge.
c::9010.20, GLM, Creation of subroutine bcorb.
c::9906.18, SAH, Convert Fortran 77 to C .
c::9906.24, SAH, added the satellite clock offset computation.
c********1*********2*********3*********4*********5*********6*********7**
ce::bcorb
c
*/
int bcorb(double tc, int isv, double recf[4], double vecf[3]) {
  double recf1[3];
  double recf2[3];
  double vecf1[3];
  double vecf2[3];
  double tcx;
  double tdif;
  double tk;
  double w1;
  double w2;
  double xcot;
  double crit;
  int ierr;
  int n;
  int n1;
  int n2;
  int save_n1;

  int i;

  /*......  implicit initialization   */

  for (n = 0; n < 3; n++) {
    recf[n] = 0.0;
    vecf[n] = 0.0;
  }
  recf[3] = 0.0;

  /*...... 1.0  check if the satellite has orbit parameters  */

  if (nxco[isv] == 0) {
    ierr = isv;
    return (ierr);
  }

  /*...... 1.1  check current time, tc, against times of orbit elements
    and select elements before and after tc */

  save_n1 = -9999;
  n1 = 0;
  for (n = 0; n < nxco[isv]; n++) {
    xcot = xco[isv][8][n];
    tk = tc - xcot;
    if (tk > 302400.0)
      xcot = xcot + 604800.0;
    else if (tk < -302400.0)
      xcot = xcot - 604800.0;

    if (xcot > tc) {
      save_n1 = n1;
      break;
    } else
      n1 = n;
  }

  if (save_n1 == -9999)
    n1 = nxco[isv] - 1;

  n2 = n1 + 1;
  if (n2 > (nxco[isv] - 1))
    n2 = nxco[isv] - 1;
  xcot = xco[isv][8][n1];
  tk = tc - xcot;
  if (tk > 302400.0)
    xcot = xcot + 604800.0;
  else if (tk < -302400.0)
    xcot = xcot - 604800.0;

  if (xcot > tc)
    n2 = n1;
  xcot = xco[isv][8][n2];
  tk = tc - xcot;
  if (tk > 302400.0)
    xcot = xcot + 604800.0;
  if (tk < -302400.0)
    xcot = xcot - 604800.0;
  if (xcot < tc)
    n1 = n2;
  if (tc == 604800.0)
    tc = 0.0;

  /*...... 2.0  compute ephemeris for this set of elements   */

  tcx = tc;
  jxco = n1;
  ierr = bccalc(tcx, isv, recf1, vecf1);
  if (ierr != 0) {
    printf("\nERROR returned from jxco=n1 bccalc for tcx = %lf isv: %2d\n", tcx,
           isv);
    return (ierr);
  }

  /*...... 3.0  compute ephemeris for this set of elements   */

  tcx = tc;
  jxco = n2;
  ierr = bccalc(tcx, isv, recf2, vecf2);
  if (ierr != 0) {
    cout << "\nERROR returned from jxco=n2 bccalc for tcx = " << tcx
         << " isv: " << isv << endl;
    return (ierr);
  }

  /*...... 4.0  compute weighted average of positions and velocities   */

  if (n1 == n2) {
    w1 = 0.50;
    w2 = 0.50;
  } else {
    tdif = xco[isv][8][n2] - xco[isv][8][n1];
    w1 = (tdif - (tc - xco[isv][8][n1])) / tdif;
    w2 = (tdif - (xco[isv][8][n2] - tc)) / tdif;
  }

  /* broadcast blocks should be used within 2 hours of TOE */
  /* record/flag every thing beyond 4 hours + 60 seconds */
  if (fabs(tc - xco[isv][8][n1]) > 14460 ||
      fabs(xco[isv][8][n2] - tc) > 14460) {
    cuml_bad_blocks[isv] = cuml_bad_blocks[isv] + 1;
  }

  for (n = 0; n < 3; n++) {
    recf[n] = w1 * recf1[n] + w2 * recf2[n];
    vecf[n] = w1 * vecf1[n] + w2 * vecf2[n];
  }

  /*...... 5.0  compute sv clock offset using the closest PRN block */

  crit = 1.0e308;
  for (n = 0; n < nxco[isv]; n++) // find the closest svclk coeffs.
  {
    xcot = xco[isv][24][n];
    if (fabs(tc - xcot) < crit) {
      crit = fabs(tc - xcot);
      tdif = tc - xcot;
      n1 = n;
    }
  }
  recf[3] = xco[isv][25][n1] + xco[isv][26][n1] * tdif +
            xco[isv][27][n1] * tdif * tdif;

  return (0);

} /* end of function bcorb */

/*
c
c********1*********2*********3*********4*********5*********6*********7**
c
c name:        bccalc
c version:     0007.19
c written by:  Steve Hilla
c purpose:     This subroutine uses the broadcast ephemeris
c              representation parameters to compute the position
c              vector (x,y,z) based on the ICD-200 (01 July 1992)
c
c input parameters
c ------------------------------------
c isv               the satellite PRN id number
c tsubr             GPS time of current epoch [seconds]
c
c output parameters
c -------------------------------------
c ierr              error flag
c                     = 1 => there are no sv's
c recf              position vector (x,y,z) and sv clock offset
c vecf              velocity vector (xdot, ydot, zdot)
c
c local variables and constants:
c -------------------------------------
c a                 semi-major axis [meters]
c clkoff            satellite clock offset
c cori              2nd harmonic correction to inclination
c corid             derivative of cori with respect to time
c corlat            2nd harmonic argument of latitude correction
c corlatd           derivative of corlat with respect to time
c corr              2nd harmonic radius correction
c corrd             derivative of corr with respect to time
c cos2p             intermediate result (cos twolat)
c cos2pd            derivative of cos2p with respect to time
c cose              cosine of eccentric anomaly (ek) [n.u.]
c cosi              cosine of xik
c cosid             derivative of cosi with respect to time
c coslat            intermediate result (cos dlat)
c coso              cosine of omegak
c dlat              argument of latitude (phi)
c                     also corrected argument of latitude
c dlatd             derivative of dlat with respect to time
c ek                eccentric anomaly
c ekdot             derivative of eccentric anomaly (ek)  w.r.t.  time
c emk               mean anomaly, emk= em0() + en*tk
c en                corrected mean motion, en= en0 + deltan() [seconds]
c en0               computed mean motion, en0= sqrt(xmu/a**3) [seconds]
c ghadot            angular velocity of the earth [radians/seconds]
c                   = 7.292115147d-5 radians/seconds
c omegak            corrected longitude of ascending node
c omegakd           derivative of omegak with respect to time
c                   (rate of right ascension)
c pi                numerical constant for pi [n.u.]
c rk                corrected radius
c rkd               derivative of rk with respect to time
c roote             intermediate result (sqrt(1-e**2))
c sin2p             intermediate result (sin twolat)
c sin2pd            derivative of sin2p with respect to time
c sine              sine of eccentric anomaly (ek) [n.u.]
c sini              sine of xik
c sinid             derivative of sini with respect to time
c sino              sine of omegak
c sinlat            hold intermediate result (sin dlat)
c truea             true anomaly
c tk                actual total time difference between the GPS system
c                   time and the ephemeris reference epoch time:
c                     tk= tc - toe()
c twolat            intermediate result (2*argument of latitude),
c                     twolat= mod(dlat+dlat/twopi)
c twopi             hold pi*pi, twopi= pi + pi, [no units]
c xik               corrected inclination
c xikd              derivative of xik with respect to time
c xk                x position in orbital plane, xk= rk*coslat
c xkd               derivative of xk with respect to time
c xmu               WGS-72 value of the earth's
c                   universal gravitational parameter
c                   = 3.986008d+14 [meters**3/seconds**2]
c yk                y position in orbital plane, yk= rk*sinlat
c ykd               derivative of yk with respect to time
c
c
c The following global variables are the ephemeris representation
c model which is characterized by a set of parameters that is an
c extension to the Keplerian orbital parameters, one set for each sv.
c
c cic               amplitude of the cosine harmonic correction term to
c                   the angle of inclination
c cis               amplitude of the sine harmonic correction term to
c                   the angle of inclination
c crc               amplitude of the cosine harmonic correction term to
c                   the orbit radius
c crs               amplitude of the sine harmonic correction term to
c                   the orbit radius
c cuc               amplitude of the cosine harmonic correction term to
c                   the argument of latitude
c cus               amplitude of the sine harmonic correction term to
c                   the argument of latitude
c deltan            delta(n) : mean motion difference from computed
c e                 value eccentricity
c em0               mean anomaly at reference time toe
c omega0            right ascension at reference time toe
c omegad            rate of right ascension
c per0              argument of perigee
c roota             square root of the semi-major axis (a)
c toe               ephemeris reference time
c xi0               inclination angle at reference time toe
c
c global variables and constants:
c -------------------------------------
c jxco              broadcast element set selected
c xco()             broadcast element sets
c
c calls:
c -------------------------------------
c
c include files:
c -------------------------------------
c
c common blocks:
c -------------------------------------
c
c references:
c -------------------------------------
c ICD-GPS-200 Interface Control Document
c IRN-200B-PR-001
c 01 July 1992
c
c comments:
c -------------------------------------
c All angular values must have been converted to radians.
c
c********1*********2*********3*********4*********5*********6*********7**
c::modification history:
c::0007.19, SAH, create bccalc using the ICD-200 and 1995 class notes.
c********1*********2*********3*********4*********5*********6*********7**
ce::bccalc
c
*/
int bccalc(double tsubr, int isv, double recf[3], double vecf[3]) {
  double a;
  double cic;
  double cis;
  double cori;
  double corid;
  double corlat;
  double corlatd;
  double corrd;
  double corr;
  double cos2p;
  double cos2pd;
  double cose;
  double cosi;
  double coslat;
  double coso;
  double crc;
  double crs;
  double cuc;
  double cus;
  double deltan;
  double dlat;
  double dlatd;
  double e;
  double ek;
  double ekdot;
  double em0;
  double emk;
  double en;
  double en0;
  double omega0;
  double omegad;
  double omegak;
  double omegakd;
  double per0;
  double rk;
  double rkd;
  double roota;
  double roote;
  double sin2p;
  double sin2pd;
  double sine;
  double sini;
  double sinlat;
  double sino;
  double tk;
  double toe;
  double truea;
  double twolat;
  double xi0;
  double xik;
  double xikd;
  double xk;
  double xkd;
  double yk;
  double ykd;
  double idot, toc, a0, a1, a2, clkoff, Ecc_anom;
  // double  relCorr;

  int i;
  int ierr;
  double xghadot = 7.2921151467e-05; // ICD-GPS-200  01 July 1992

  /*....... 1.0  check to see if isv exists  */

  // cout << " isv,jxco,tc = " << isv << " " << jxco << " " <<
  //  setw(30) << setprecision(15) << tsubr << endl;

  if (isv == 0) {
    ierr = 1;
    return (ierr);
  }

  /*....... 2.0  compute the earth-centered,earth-fixed cartesian
    coordinates of the position for isv's antenna phase
    center  */

  crs = xco[isv][1][jxco];
  deltan = xco[isv][2][jxco];
  em0 = xco[isv][3][jxco];

  cuc = xco[isv][4][jxco];
  e = xco[isv][5][jxco];
  cus = xco[isv][6][jxco];
  roota = xco[isv][7][jxco];

  toe = xco[isv][8][jxco];
  cic = xco[isv][9][jxco];
  omega0 = xco[isv][10][jxco];
  cis = xco[isv][11][jxco];

  xi0 = xco[isv][12][jxco];
  crc = xco[isv][13][jxco];
  per0 = xco[isv][14][jxco];
  omegad = xco[isv][15][jxco];

  idot = xco[isv][16][jxco];

  //      toc=    xco[isv][24][jxco];
  //      a0=     xco[isv][25][jxco];
  //      a1=     xco[isv][26][jxco];
  //      a2=     xco[isv][27][jxco];

  // cout << " toe = " << setw(30) << setprecision(15) << toe << endl;
  // cout << " omega0 = " << setw(30) << setprecision(15) << omega0 << endl;
  // cout << " per0 = " << setw(30) << setprecision(15) << per0 << endl;
  // cout << " omegad = " << setw(30) << setprecision(15) << omegad << endl;
  // cout << " xi0 = " << setw(30) << setprecision(15) << xi0 << endl;
  // cout << " em0 = " << setw(30) << setprecision(15) << em0 << endl;

  /*....... 2.01  compute SV clock offset and relativity correction
    apply these to the input transmission time */

  //      if( (tsubr - toc) > 302400.0 ) tsubr = tsubr - 604800.0;
  //      if( (tsubr - toc) < -302400.0 ) tsubr = tsubr + 604800.0;
  //      clkoff = a0 + a1*(tsubr - toc) + a2*(tsubr - toc)*(tsubr - toc);

  //      Ecc_anom = em0;
  //      for( i = 0; i < 7; i++ )
  //        Ecc_anom = em0 + e*sin(Ecc_anom);

  //      relCorr = -4.442807633e-10 * e * roota * sin(Ecc_anom);
  //
  //      tsubr = tsubr - clkoff - relCorr;

  /*....... 2.1  compute semi-major axis and computed mean motion  */

  a = roota * roota;
  en0 = sqrt(xmu / (a * a * a));

  /*....... 2.2  compute time from epoch and account for beginning or
    end of week crossovers.  */

  tk = tsubr - toe;
  if (tk > 302400.0)
    tk = tk - 604800.0;
  else if (tk < -302400.0)
    tk = tk + 604800.0;

  /*....... 2.3  compute corrected mean motion and mean anomaly  */

  en = en0 + deltan;
  emk = em0 + en * tk;
  emk = fmod(emk, twopi);
  emk = fmod(emk + twopi, twopi);

  /*....... 2.4  solve Kepler's equation for eccentric anomaly.

    note that Kepler's equation   emk= ek-e*sin(ek)
    can be re-written in the form  ek= emk+e*sin(ek)

    ek= emk is the initial approximation for Kepler's
    equation and ekdot is the derivative of ek with
    respect to time.   */

  ek = emk;
  ekdot = en;

  for (i = 0; i < 4; i++) {
    sine = sin(ek);
    cose = cos(ek);
    ek = emk + e * sine;
    ekdot = en + e * cose * ekdot;
  }

  /*....... 2.5  compute true anomaly, which is derived the following way:

    sin(vk)   sqrt(1-e**2)*sin(ek)   1-e*cos(ek)
    tan(vk)= -------= -------------------- * -----------
    cos(vk)       1-e*cos(ek)         cos(ek)-e  */

  if ((1.0 - e * e) < 1.0e-4) {
    cout << " 1.0 - e*e is less than 0.0001 ! in bccalc.\n";
  }
  roote = sqrt(1.0 - e * e);
  truea = atan2(roote * sine, cose - e);

  /*....... 2.6  compute argument of latitude, and 2nd harmonic pertrbtns.
    which include argument of latitude correction (corlat),
    radius correction (corr), and corr. to inclination (cori)  */

  dlat = truea + per0;
  twolat = fmod(dlat + dlat, twopi);
  sin2p = sin(twolat);
  cos2p = cos(twolat);

  corlat = cuc * cos2p + cus * sin2p;
  corr = crc * cos2p + crs * sin2p;
  cori = cic * cos2p + cis * sin2p;

  /*....... 2.7  compute corrected argument of latitude, radius, and
    inclination (dlat, rk, and xik) respectively  */

  dlat = fmod(dlat + corlat, twopi);
  rk = a * (1.0 - e * cose) + corr;
  xik = xi0 + idot * tk + cori;

  /*....... 2.8  compute x and y positions in orbital plane  */

  coslat = cos(dlat);
  sinlat = sin(dlat);
  xk = rk * coslat;
  yk = rk * sinlat;
  // cout << " rk = " << setw(30) << setprecision(15) << rk << endl;
  // cout << " dlat = " << setw(30) << setprecision(15) << dlat << endl;

  /*....... 2.9  compute corrected longitude of ascending node (omegak),
    and finally, compute the earth-centered, earth-fixed
    cartesian coordinates of the position for isv's antenna
    phase center (recf)  */

  omegak = omega0 + (omegad - xghadot) * tk - (xghadot * toe);
  omegak = fmod(omegak, twopi);
  coso = cos(omegak);
  sino = sin(omegak);
  cosi = cos(xik);
  sini = sin(xik);
  recf[0] = xk * coso - yk * cosi * sino;
  recf[1] = xk * sino + yk * cosi * coso;
  recf[2] = yk * sini;

  //      double tdif = tsubr - xco[isv][24][jxco];
  //     recf[3] = xco[isv][25][jxco] + xco[isv][26][jxco]*tdif
  //                + xco[isv][27][jxco]*tdif*tdif;

  // cout << " xghadot = " << setw(30) << setprecision(15) << xghadot << endl;
  // cout << " xk = " << setw(30) << setprecision(15) << xk << endl;
  // cout << " yk = " << setw(30) << setprecision(15) << yk << endl;
  // cout << " coso = " << setw(30) << setprecision(15) << coso << endl;
  // cout << " sino = " << setw(30) << setprecision(15) << sino << endl;
  // cout << " cosi = " << setw(30) << setprecision(15) << cosi << endl;
  // cout << " sini = " << setw(30) << setprecision(15) << sini << endl;
  // cout << " x = " << setw(30) << setprecision(15) << recf[0] << endl;
  // cout << " y = " << setw(30) << setprecision(15) << recf[1] << endl;
  // cout << " z = " << setw(30) << setprecision(15) << recf[2] << endl;
  //  cin >> i;

  /*....... 4.1  compute the derivative of omegak and dlat w.r.t. time tk.
    note:  the formula for the derivative of dlat is
    fairly complicated to work out.  */

  omegakd = omegad - ghadot;
  dlatd = en * roote * ((a / rk) * (a / rk));

  /*....... 4.2  compute the derivative of the 2nd harmonic corrections
    with respect to time tk. (i.e. compute the derivative of:
    1 - argument of latitude correction
    2 - radius correction
    3 - correction to inclination)  */

  sin2pd = 2.0 * cos2p * dlatd;
  cos2pd = -2.0 * sin2p * dlatd;
  corlatd = cuc * cos2pd + cus * sin2pd;
  corrd = crc * cos2pd + crs * sin2pd;
  corid = cic * cos2pd + cis * sin2pd;

  /*....... 4.3  compute the derivative of the following w.r.t. time tk:
    1) corrected argument of latitude
    2) corrected radius
    3) corrected inclination  */

  dlatd = dlatd + corlatd;
  rkd = a * e * sine * ekdot + corrd;
  xikd = corid;

  /*....... 4.4  compute the derivative of the x and y positions in the
    orbital plane, and finally, compute the velocity vector,
    which is the derivative of the earth's fixed coordinates  */

  xkd = rkd * coslat - rk * sinlat * dlatd;
  ykd = rkd * sinlat + rk * coslat * dlatd;

  vecf[0] = (xkd * coso) - (ykd * cosi * sino) - (xk * sino * omegakd) -
            (yk * cosi * coso * omegakd) + (yk * sini * xikd * sino);

  vecf[1] = (xkd * sino) + (ykd * cosi * coso) + (xk * coso * omegakd) -
            (yk * cosi * sino * omegakd) - (yk * sini * xikd * coso);

  vecf[2] = (ykd * sini) + (yk * cosi * xikd);

  ierr = 0;
  return (ierr);

} /* end of function bccalc */

//---------------------------------------------------------------------------

int bcread(RinexNavFile &navFile, ofstream &out) {
  double tel;
  int n, isv;
  int duplicate_block;

  PRNBlock currBlock;
  DateTime currEpoch;
  MJD modJD;
  GPSTime gpsTime;

  /*.......  1.0  explicit initialization */

  for (isv = 0; isv < MAXSVS; isv++) {
    nxco[isv] = 0;
  }
  /*.......  2.0  read the broadcast header */

  navFile.readHeader();
  navFile.writeHeaderImage(out);
  out << endl << endl;

  /*.......       TOP OF LOOP to scan file   */
  /*.......  3.0  read epoch / clock model line   */

  while (navFile.readPRNBlock(currBlock) != 0) {
    //     navFile.writePRNBlock(out, currBlock);
    isv = currBlock.getSatellitePRN();

    /*.......  4.0  convert time to DateTime   */

    currEpoch =
        DateTime((long)currBlock.getTocYear(), (long)currBlock.getTocMonth(),
                 (long)currBlock.getTocDay(), (long)currBlock.getTocHour(),
                 (long)currBlock.getTocMin(), currBlock.getTocSec());
    modJD = currEpoch.GetMJD();
    tel = (double)modJD.mjd + modJD.fracOfDay;

    /*.......  5.0  compute GPS Time   */

    gpsTime = currEpoch.GetGPSTime();

    /*.......  6.0  skip blocks that are out of order, go to end of while loop
     */

    if (nxco[isv] >= 1 && tel < xco[isv][20][nxco[isv] - 1])
      continue;

    /*.......  6.1  skip blocks where the PRN is declared unhealthy
      (zero means the Satellite Health is OK) */

    if (fabs(currBlock.getSvHealth()) > 0.00000000000001) {
      out << "In the NAV file, PRN " << setw(2) << isv
          << " is unhealthy at time: " << setw(4) << currBlock.getTocYear()
          << " " << setw(2) << currBlock.getTocMonth() << " " << setw(2)
          << currBlock.getTocDay() << " " << setw(2) << currBlock.getTocHour()
          << " " << setw(2) << currBlock.getTocMin() << fixed << showpoint
          << setw(7) << setprecision(3) << currBlock.getTocSec() << endl;
      continue;
    }

    /*.......  7.0  skip duplicate blocks
      NOTE: only the epoch is checked.  */

    duplicate_block = 0;
    if (nxco[isv] > 1) {
      for (n = 0; n < nxco[isv] - 2; n++) {
        if (fabs(xco[isv][8][n] - currBlock.getToe()) < 1.0) {
          duplicate_block = 1;
          break; /* exit the for loop */
        }
      }

      if (duplicate_block == 1)
        continue;
    }

    /*.......  8.0  load elements into the storage array  */

    if (nxco[isv] < MAXEPOCH) {
      xco[isv][0][nxco[isv]] = currBlock.getIode();
      xco[isv][1][nxco[isv]] = currBlock.getCrs();
      xco[isv][2][nxco[isv]] = currBlock.getDeltan();
      xco[isv][3][nxco[isv]] = currBlock.getMo();

      xco[isv][4][nxco[isv]] = currBlock.getCuc();
      xco[isv][5][nxco[isv]] = currBlock.getEccen();
      xco[isv][6][nxco[isv]] = currBlock.getCus();
      xco[isv][7][nxco[isv]] = currBlock.getSqrtA();

      xco[isv][8][nxco[isv]] = currBlock.getToe();
      xco[isv][9][nxco[isv]] = currBlock.getCic();
      xco[isv][10][nxco[isv]] = currBlock.getBigOmega();
      xco[isv][11][nxco[isv]] = currBlock.getCis();

      xco[isv][12][nxco[isv]] = currBlock.getIo();
      xco[isv][13][nxco[isv]] = currBlock.getCrc();
      xco[isv][14][nxco[isv]] = currBlock.getLilOmega();
      xco[isv][15][nxco[isv]] = currBlock.getBigOmegaDot();

      xco[isv][16][nxco[isv]] = currBlock.getIdot();
      xco[isv][17][nxco[isv]] = currBlock.getCodesOnL2();
      xco[isv][18][nxco[isv]] = currBlock.getToeGPSWeek();
      xco[isv][19][nxco[isv]] = currBlock.getPDataFlagL2();

      xco[isv][20][nxco[isv]] = currBlock.getSvAccur();
      xco[isv][21][nxco[isv]] = currBlock.getSvHealth();
      xco[isv][22][nxco[isv]] = currBlock.getTgd();
      xco[isv][23][nxco[isv]] = currBlock.getIodc();

      xco[isv][24][nxco[isv]] = gpsTime.secsOfWeek;
      xco[isv][25][nxco[isv]] = currBlock.getClockBias();
      xco[isv][26][nxco[isv]] = currBlock.getClockDrift();
      xco[isv][27][nxco[isv]] = currBlock.getClockDriftRate();

      xco[isv][20][nxco[isv]] = tel; /* use element 20 to store mjd + fmjd */
      nxco[isv] = nxco[isv] + 1;
    }

  } /* end of while loop over all the broadcast blocks */

  return (0);

} /* end of bcread */

// -------------------------------------------------------------------

string padZeros(int a) {
  ostringstream stream;
  string tmp;
  stream << a;
  // cout <<  "length in padZeros " << stream.str().length() << endl;
  if (stream.str().length() == 1) {
    tmp = "0" + stream.str();
  } else {
    tmp = stream.str();
  }

  return tmp;
}

/* ---------------------------------------------- */
/* Anderson & Mikhail, Surveying Theory & Practice */
/* pg 1122 */

void llh2xyz(double a, double finverse, double lat, double lon, double h,
             double *x, double *y, double *z) {
  int iter;
  double e2, delta, tmp, N, lat1, lat_prev, htmp, f;
  f = 1.0 / finverse;
  e2 = (2.0 - f) * f;

  N = a / sqrt(1.0 - e2 * sin(lat) * sin(lat));

  *x = (N + h) * cos(lat) * cos(lon);
  *y = (N + h) * cos(lat) * sin(lon);
  *z = (N * (1.0 - e2) + h) * sin(lat);
}
