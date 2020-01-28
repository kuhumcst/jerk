/* JERK
Analyses the movement of two points in x-y plane, in casu nose tips data from
OpenPoseDemo.exe, and computes velocity, acceleration and jerk of the points.

The program expects input with data for two persons in a dialogue situation.

The format of the input is as follows:

FrameNr <tab> <person A data> <tab> <person B data> <newline>
  ...

where data for each person is

<x-position> <tab> <y-position> <tab> <probability> 

The probability value is currently ignored.

The command line arguments to the program are
  input file
  output file
  "name" of the first person
  "name" of the second person
  number of observations used for computing velocity
  number of observations used for computing acceleration
  number of observations used for computing jerk
  location of the observations with respect to the created annotation: 
       'past' or 'middle'
  
For example

  jerk F2_M4_SPLIT_FINAL.tab F2_M4_SPLIT_FINAL-7-14-21-past.vaj.tab F2 M4 7 14 21 past

  
Output format is
<frame number x 4> <tab> <person A data> <tab> <person B data> <new line>
  ...

<person data> = 
<x> <tab> <y> <tab> <probability> <tab> <velocity> <tab> <acceleration> <tab> <jerk> <new line>

where
  <velocity> = <strength> <tab> <angle (1..12)> <tab> <x component> <tab> <y component>
   etc.

<strength> and <angle> are just another (polar) representation of the components.
1..12 are 'clock' directions. '12' is up, '3' to the right, etc.

Author: Bart Jongejan
2020.01.23
*/

/*  */

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define M_PI           3.14159265358979323846  /* pi */

enum
    {
    evelo, eacce, ejerk
    };

struct person
    {
    double x;
    double y;
    double p;
    double vx;
    double ax;
    double jx;
    double vy;
    double ay;
    double jy;
    };

struct line
    {
    int frame;
    struct person P[2];
    };

double velocity(double St, double St2, double Sh, double Sth, double period)
    {
/*
first degree: a:position, b:velocity
var a
solution (a.period^-1*(Sh+-1*St*b))
var b
  solution
  (b.(-1*St^2+St2*period)^-1*(-1*Sh*St+Sth*period))

For vertical, replace h th by v tv
*/
    return (Sth * period - Sh * St) / (St2 * period - St * St);
    }

double acceleration(double St, double St2, double St3, double St4, double Sh, double Sth, double St2h, double period)
    {

/*

acceleration c

(horizontal)

accumulate (S) t t2 t3 t4 h th t2h
(vertical                 v tv t2v )

solution
( c
.     ( -1*St2^3
    + 2*St*St2*St3
    + St2*St4*period
    + -1*St^2*St4
    + -1*St3^2*period
    )
  ^ -1
* ( -1*Sh*St2^2
  + Sh*St*St3
  + St*St2*Sth
  + St2*St2h*period
  + -1*St3*Sth*period
  + -1*St^2*St2h
  )
)

For vertical, replace h th t2h by v tv t2v
*/
    return
        (Sh * (St * St3 - St2 * St2)
         + Sth * (St * St2 - St3 * period)
         + St2h * (St2 * period - St * St)
         )
        /
        (St2 * (2 * St * St3 - St2 * St2 + St4 * period)
         - St * St * St4
         - St3 * St3 * period
         );

    }

double jerk(double St2, double St2h, double St3, double St3h, double St4, double St5, double St6, double Sth, double period)
    {
    double var1 = St3 * St3;
    double var2 = St4 * St4;
    double var5 = St3 * St5;
    double var6 = var1 - St2 * St4;

    return -(St2h * (St2 * St5 * period - St3 * (St2 * St2 + St4 * period)) + St3h * (St2 * St2 * St2 + period * var6) + Sth * (St2 * var6 + period * (var2 - var5)))
        / (var1 * var1 - St2 * (St5 * St5 * period - St2 * (var2 + 2 * var5 - St2 * St6) + St4 * (3 * var1 - St6 * period)) - period * (St6 * var1 + St4 * (var2 + -2.0 * var5)));
    }


void doPerson(int Person, int what, const char* Period, const char* margin, line* Lines)
    {
    size_t seqsiz = strtoul(Period, 0, 10); 
    double period = (double)seqsiz;
    long t0;
    long* Ts = (long*)calloc(seqsiz, sizeof(long));
    double* ts = (double*)calloc(seqsiz, sizeof(double));
    double* hs = (double*)calloc(seqsiz, sizeof(double));
    double* vs = (double*)calloc(seqsiz, sizeof(double));
    size_t index = 0;
    struct line* current = Lines;

    long Offset;
    if(!strcmp(margin, "middle"))
        Offset = (seqsiz -1) >> 1;
    else
        Offset = 0;
    for(current = Lines; current->frame != 0; ++current)
        {
        struct person* pers = current->P + Person;
        struct line* assignTo = current - Offset;
        struct person* assignToPers = assignTo->P + Person;
        if(assignTo < Lines)
            assignTo = 0;
        else
            switch(what)
                {
                case evelo:
                    assignToPers->vx = 0;
                    assignToPers->vy = 0;
                    break;
                case eacce:
                    assignToPers->ax = 0;
                    assignToPers->ay = 0;
                    break;
                case ejerk:
                    assignToPers->jx = 0;
                    assignToPers->jy = 0;
                    break;
                }

        if(pers->x != 0)
            {
            int ind = index % seqsiz;

            if(index == 0)
                // Start accumulation of datapoints
                {
                t0 = current->frame;
                Ts[ind] = 0;
                }
            else
                {
                Ts[ind] = current->frame - t0;
                }
            hs[ind] = pers->x;
            vs[ind] = pers->y;
            ++index;
            if(index >= seqsiz)
                {
                size_t i;
                double averageTime = 0;

                double St = 0.0, St2 = 0.0, St3 = 0.0, St4 = 0.0, St5 = 0.0, St6 = 0.0,
                    Sh = 0.0, Sth = 0.0, St2h = 0.0, St3h = 0.0,
                    Sv = 0.0, Stv = 0.0, St2v = 0.0, St3v = 0.0;
                double averageh = 0;
                double averagev = 0;
                long Tbias = 0; 
                for(i = 0; i < seqsiz; ++i)
                    {
                    Tbias += Ts[i];
                    }
                Tbias /= seqsiz;
                for(i = index - seqsiz; i < index; ++i)
                    {
                    int id = i % seqsiz;
                    averageh += hs[id];
                    averagev += vs[id];
                    ts[id] = (double)(Ts[id]-Tbias); // This makes ts[id] a small number around 0
                    averageTime += ts[id]; // later used to slightly adjust the ts[id] so that their sum becomes zero.
                    }
                averageh = averageh / period;
                averagev = averagev / period;
                averageTime = averageTime / period;
                for(i = index - seqsiz; i < index; ++i)
                    {
                    int id = i % seqsiz;
                    double t = ts[id] - averageTime;
                    double t2 = t * t;
                    double t3 = t2 * t;
                    double h = hs[id] - averageh;
                    double v = vs[id] - averagev;
                    St += t;
                    St2 += t2;
                    St3 += t2 * t;
                    St4 += t2 * t2;
                    St5 += t2 * t3;
                    St6 += t3 * t3;
                    Sh += h;
                    Sth += t * h;
                    St2h += t2 * h;
                    St3h += t3 * h;
                    Sv += v;
                    Stv += t * v;
                    St2v += t2 * v;
                    St3v += t3 * v;
                    }
                switch(what)
                    {
                    case evelo:
                        assignToPers->vx = velocity(St, St2, Sh, Sth, period);
                        assignToPers->vy = velocity(St, St2, Sv, Stv, period);
                        break;
                    case eacce:
                        assignToPers->ax = acceleration(St, St2, St3, St4, Sh, Sth, St2h, period);
                        assignToPers->ay = acceleration(St, St2, St3, St4, Sv, Stv, St2v, period);
                        break;
                    case ejerk:
                        assignToPers->jx = jerk(St2, St2h, St3, St3h, St4, St5, St6, Sth, period);
                        assignToPers->jy = jerk(St2, St2v, St3, St3v, St4, St5, St6, Stv, period);
                        break;
                    }
                }
            else
                {
                //Not enough data collected to fill a period.
                }
            }
        }
    }

line* readInput(char* name)
    {
    FILE* fp;
    printf("name %s\n", name);
    fp = fopen(name, "r");
    if(fp)
        {
        int lines = 1;
        int kar;
        while((kar = fgetc(fp)) != EOF)
            {
            if(kar == '\n')
                ++lines;
            }
        rewind(fp);
        // 140	94.7639	90.7244	0.357086	590.125	110.799	0.275333
        line* Lines = (line*)calloc(lines, sizeof(struct line));
        int L = 0;
        int frame;
        double xl;
        double yl;
        double pl;
        double xr;
        double yr;
        double pr;
        for(L = 0; L < lines -1 && fscanf(fp, "%d\t%lg\t%le\t%lf\t%lf\t%lf\t%lf\r\n", &frame, &xl, &yl, &pl, &xr, &yr, &pr); ++L)
            {
            Lines[L].frame = frame;
            Lines[L].P[0].x = xl;
            Lines[L].P[0].y = yl;
            Lines[L].P[0].p = pl;
            Lines[L].P[1].x = xr;
            Lines[L].P[1].y = yr;
            Lines[L].P[1].p = pr;
            }
        Lines[L].frame = 0;
        Lines[L].P[0].x = 0.0;
        Lines[L].P[0].y = 0.0;
        Lines[L].P[0].p = 0.0;
        Lines[L].P[1].x = 0.0;
        Lines[L].P[1].y = 0.0;
        Lines[L].P[1].p = 0.0;
        fclose(fp);
        return Lines;
        }
    else
        return 0;
    }

int clock(double sagittaH, double sagittaV)
    {
    double theta = atan2((double)-sagittaV, (double)sagittaH) * 180.0 / M_PI;
    theta = -theta;
    if(theta < 0.0)
        theta += 360;
    theta += 105; // half past eleven
    int direction = (int)theta;
    direction %= 360;
    direction /= 30; // 0 .. 11
    if(direction == 0)
        direction = 12;
    return direction;
    }

double rho(double sagittaH, double sagittaV)
    {
    double rho = sqrt((double)(sagittaV * sagittaV) + (double)(sagittaH * sagittaH));
    return rho;
    }

void printPerson(FILE* fp, struct person* Person, char* name)
    {
    fprintf(fp,
            "%lf\t%lf\t%lf\t%lf\t%ld\t%lf\t%lf\t%lf\t%ld\t%lf\t%lf\t%lf\t%ld\t%lf\t%lf",
            Person->x, Person->y, Person->p,
            rho(Person->vx, Person->vy), clock(Person->vx, Person->vy), Person->vx, Person->vy,
            rho(Person->ax, Person->ay), clock(Person->ax, Person->ay), Person->ax, Person->ay,
            rho(Person->jx, Person->jy), clock(Person->jx, Person->jy), Person->jx, Person->jy
    );
            // F2:x-pos	F2:y-pos	F2:weight	F2:velocity-r	F2:velocity-clock	F2:velocity-x	F2:velocity-y	F2:acceleration-r	F2:acceleration-clock	F2:acceleration-x	F2:acceleration-y	F2:jerk-r	F2:jerk-clock	F2:jerk-x	F2:jerk-y
    }

void printPersonHead(FILE* fp, char* name)
    {
    fprintf(fp
            , "%s:x-pos\t%s:y-pos\t%s:weight\t%s:velocity-r\t%s:velocity-clock\t%s:velocity-x\t%s:velocity-y\t%s:acceleration-r\t%s:acceleration-clock\t%s:acceleration-x\t%s:acceleration-y\t%s:jerk-r\t%s:jerk-clock\t%s:jerk-x\t%s:jerk-y"
            , name, name, name, name, name, name, name, name, name, name, name, name, name, name, name
    );
    }

static char names[2][3];

void print(struct line* Lines, const char* out)
    {
    FILE* fp = fopen(out, "wb");
    if(fp)
        {
        fprintf(fp, "frame\t");
        printPersonHead(fp, names[0]);
        fputc('\t', fp);
        printPersonHead(fp, names[1]);
        fputc('\n', fp);
        struct line* current = Lines;
        for(current = Lines; current->frame; ++current)
            {
            fprintf(fp, "%d\t", (current->frame) << 2);
            for(int p = 0; p < 2; ++p)
                {
                printPerson(fp, current->P + p, names[p]);
                if(p == 0)
                    fputc('\t', fp);
                else
                    fputc('\n', fp);
                }
            }
        fclose(fp);
        }
    }

int main(int argc, char** argv)
    {
    if(argc != 9)
        {
        //       0      1       2          3             4               5                   6                 7             8
        printf("jerk <input> <output> <leftperson> <rightperson> <velocity period> <acceleration period> <jerk period> <margin type>\n");
        return -1;
        }
    names[0][0] = argv[1][0];
    names[0][1] = argv[1][1];
    names[0][2] = 0;
    names[1][0] = argv[1][3];
    names[1][1] = argv[1][4];
    names[1][2] = 0;
    struct line* Lines = readInput(argv[1]);
    if(Lines)
        {
        unsigned long velow = strtoul(argv[5], NULL, 10);
        unsigned long accew = strtoul(argv[6], NULL, 10);
        unsigned long jerkw = strtoul(argv[7], NULL, 10);
        int pers = 0;
        for(pers = 0; pers < 2; ++pers)
            {
            doPerson(pers, evelo, argv[5], argv[8], Lines);
            doPerson(pers, eacce, argv[6], argv[8], Lines);
            doPerson(pers, ejerk, argv[7], argv[8], Lines);
            }
        print(Lines, argv[2]);
        free(Lines);
        }
    return 0;
    }
