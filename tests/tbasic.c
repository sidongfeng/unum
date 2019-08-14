/*
$Id: tbasic.c $

Description: Test file

$Log: $
*/

#include <stdio.h> /* printf, putchar */
#include <stdlib.h> /* exit, atoi, malloc, realloc, free, rand */
#include <string.h> /* strchr, strlen, strcmp, memcpy */
#include <ctype.h> /* isdigit */
#include <errno.h> /* errno */
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
//#define NDEBUG 1
//#include <assert.h> /* assert */

#include "unum.h"
#include "../src/conv.c"

#define PROG "tbasic"
#ifndef VERSION
#define VERSION "1.0"
#endif

#define DEFAULT_INT 1
#define DEFAULT_STR "ABC"

#define BFLAG 0x01

#define VFLAG 0x1000
const unsigned int ITERS = 10000;

int flags; /* argument flags */
int iarg = DEFAULT_INT; /* int argument */
char *sarg = DEFAULT_STR; /* string argument */
int64_t TimeInMicros() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
}

int num_same(const char *a, const char *ae, const char *b, const char *be, size_t length)
{
	/* mantissa */
	while (a < ae && b < be && length) {
		if (*a != *b) return 0;
		if (isdigit(*a)) length--;
		a++; b++;
	}
	/* exponent */
	if ((a = strpbrk(a, "Ee")) != NULL && a < ae) {
		if ((b = strpbrk(b, "Ee")) != NULL && b < be) {
			a++; b++;
			while (a < ae && b < be) {
				if (*a != *b) return 0;
				a++; b++;
			}
		} else return 0;
	}
	return 1;
}

int ustr_same(const char *a, const char *b, size_t length)
{
	char *aL = strpbrk(a, "([");
	char *aC = strchr(a, ',');
	char *aR = strpbrk(a, ")]");
	char *bL = strpbrk(b, "([");
	char *bC = strchr(b, ',');
	char *bR = strpbrk(b, ")]");
	if (aL != NULL && aC != NULL && aR != NULL) { /* bound */
		if (bL != NULL && bC != NULL && bR != NULL) {
			if (*aL == *bL && *aC == *bC && *aR == *bR) {
				return
					num_same(aL+1, aC, bL+1, bC, length) &&
					num_same(aC+1, aR, bC+1, bR, length);
			}
		}
	} else { /* exact */
		return
			num_same(a, strchr(a, '\0'), b, strchr(a, '\0'), length);
	}
	return 0;
}

/* arbitrary digits of integer */
void randomInt(char randomString[],size_t length) {
    static char charset[] = "0123456789";
	
    if (length) {
        if (randomString) {
            int l = (int) (sizeof(charset) -1);
            for (int n = 0;n < length;n++) {
                int key = rand() % l; 
                randomString[n] = charset[key];
            }
            
        }
    }
    randomString[length] = '\0';
}

/* arbitrary digits of float */
void randomFloat(char *randomString,size_t real,size_t decimal) { 
    static char charset[] = "0123456789";
    if (real) {
        if (randomString) {
            int l = (int) (sizeof(charset) -1);
            for (int n = 0;n < real;n++) {
                int key = rand() % l; 
                randomString[n] = charset[key];
            }
        }
    }
    randomString[real] = '.';
    if (decimal) {
        if (randomString) {
            int l = (int) (sizeof(charset) -1);
            for (int n = 0;n < decimal;n++) {
                int key = rand() % l;
                randomString[real+n+1] = charset[key];
            }
        }
    }
    randomString[real+decimal+1] = '\0';
}

/* float precision */
float random_float(float a){
  	return ((float)rand()/(float)(RAND_MAX)) * a;
}

/* generate a random interval */
void random_interval(char *result, float range){
	char xstr[64], ystr[64], close[1];
	int sign, l, r;
	float x, y;
	x = random_float(range);
	sign = rand() % 5;
	if (!sign) x = -x;
	while (1){
		y = random_float(range);
		sign = rand() % 5;
		if (!sign) y = -y;
		if (y > x) break;
	}
	l = rand() % 2;
	r = rand() % 2;
	if (l){
		if (r){sprintf(result,"[%f,%f]", x,y);
		}else{sprintf(result,"[%f,%f)", x,y);}
	}else{
		if (r){sprintf(result,"(%f,%f]", x,y);
		}else{sprintf(result,"(%f,%f)", x,y);}
	}
}

int main(int argc, char *argv[])
{
	int nok = 0;
	char *s;
	int tfail = 0;
	srand ( time(NULL) );

	while (--argc > 0 && (*++argv)[0] == '-')
		for (s = argv[0]+1; *s; s++)
			switch (*s) {
			case 'b':
				flags |= BFLAG;
				break;
			case 'i':
				if (isdigit((int)s[1])) iarg = atoi(s+1);
				else nok = 1;
				s += strlen(s+1);
				break;
			case 's':
				sarg = s+1;
				s += strlen(s+1);
				break;
			case 'v':
				flags |= VFLAG;
				break;
			default:
				nok = 1;
				fprintf(stderr, " -- not an option: %c\n", *s);
				break;
			}

	if (flags & VFLAG) fprintf(stderr, "%s %s\n", PROG, VERSION);
	if (nok /*|| argc < 1*/ || (argc > 0 && *argv[0] == '?')) {
		fprintf(stderr, "Usage: %s -bv -i<int> -s<str> <in_file> [<out_file>]\n", PROG);
		fprintf(stderr, "  -b  boolean argument\n");
		fprintf(stderr, "  -i  integer argument, default: %d\n", DEFAULT_INT);
		fprintf(stderr, "  -s  string argument, default: %s\n", DEFAULT_STR);
		fprintf(stderr, "  -v  version\n");
		exit(EXIT_FAILURE);
	}

#if 0
	{
		FILE *fin, *fout;

		if ((fin = fopen(argv[0], "r")) == NULL) {
			fprintf(stderr, " -- can't open file: %s\n", argv[0]);
			exit(EXIT_FAILURE);
		}
		if (argc < 2) {
			fout = stdout;
		} else if ((fout = fopen(argv[1], "w")) == NULL) {
			fprintf(stderr, " -- can't open file: %s\n", argv[1]);
			exit(EXIT_FAILURE);
		}

		/* do something */

		fclose(fin);
		fclose(fout);
	}
#endif

/*----------------------*/
/*-------- unum --------*/
/*----------------------*/

#if 1
	unum_set_env(3, 4);
	{
		UNUM_VAR(u1);
		UNUM_VAR(u2);
		UNUM_VAR(ur);
		UNUM_VAR(uc);
		int ess, fss;
		int ir, ok, fail = 0;

		unum_get_env(&ess, &fss);
		printf("\n# test unum operations, env:%d,%d #\n", ess, fss);

#define UNUM_OP(oper,op1,chk) \
		unum_set_str(u1, op1); \
		unum_##oper(ur, u1); \
		unum_set_str(uc, chk); \
		fail |= !(ok = unum_same(ur, uc)); \
		printf("%s ", ok ? "OK  " : "FAIL"); \
		printf("%s ", #oper); unum_print(u1); \
		printf(" = "); unum_print(ur); putchar('\n')
		UNUM_OP(log, "1", "0");
		UNUM_OP(log, "2", "0.6931471805");
		UNUM_OP(log, "0.23487", "-1.448723109394");
		UNUM_OP(log, "-2", "NaN");
		UNUM_OP(log, "Inf", "Inf");
		UNUM_OP(log, "0", "-Inf");
		UNUM_OP(log, "123456789", "18.63140176616801803319");
		UNUM_OP(log, "(8,9)", "(2.0794415416798359,2.1972245773362193827)");
		UNUM_OP(exp, "1", "2.718281828459045235360");
		UNUM_OP(exp, "0.23487", "1.2647443412809724235");
		UNUM_OP(exp, "-2", "0.135335283236612");
		UNUM_OP(exp, "Inf", "Inf");
		UNUM_OP(exp, "-Inf", "0");
		tfail |= fail;
	}
#endif

#if 1
	unum_set_env(3,4);
	{
		UNUM_VAR(u1);
		UNUM_VAR(u2);
		UNUM_VAR(ur);
		UNUM_VAR(uc);
		int ess, fss;
		int ir, ok, fail = 0;
		char str[64];

		unum_get_env(&ess, &fss);
		printf("\n# test unum operations, env:%d,%d #\n", ess, fss);

#define UNUM_AOP(op1,oper,op2,chk) \
		unum_set_str(u1, op1); \
		unum_set_str(u2, op2); \
		unum_##oper(ur, u1, u2); \
		unum_set_str(uc, chk); \
		fail |= !(ok = unum_same(ur, uc)); \
		printf("%s ", ok ? "OK  " : "FAIL"); \
		unum_print(u1); printf(" %s ", #oper); unum_print(u2); \
		printf(" = "); unum_print(ur); putchar('\n')
		UNUM_AOP("3",pow,"4","81");
		UNUM_AOP("2",pow,"0.5","1.41421356237309504880168");
		UNUM_AOP("-5",pow,"0.231","NaN");
		UNUM_AOP("(0,1]",pow,"2.326","(0,1]");
		UNUM_AOP("-12",pow,"2","144");
		UNUM_AOP("-Inf",pow,"2.123","NaN");

		tfail |= fail;
	}
#endif

/*----------------------*/
/*-------- ubnd --------*/
/*----------------------*/

#if 1
	unum_set_env(3, 4);
	{
		UBND_VAR(ub1);
		UBND_VAR(ub2);
		UBND_VAR(ubr);
		UBND_VAR(ubc);
		int ess, fss;
		int ir, ok, fail = 0;
		char xstr[64],ystr[64];
		int sign, xsign, ysign;
		float max = 200000.0;

		unum_get_env(&ess, &fss);
		printf("\n# test ubnd operations, env:%d,%d #\n", ess, fss);

#define UBND_OP(oper,op1,chk) \
		ubnd_set_str(ub1, op1); \
		ubnd_##oper(ubr, ub1); \
		ubnd_set_str(ubc, chk); \
		fail |= !(ok = ubnd_seq(ubr, ubc)); \
		printf("%s ", ok ? "OK  " : "FAIL"); \
		printf("%s ", #oper); ubnd_print(ub1); \
		printf(" = "); ubnd_print(ubr); putchar('\n')

		UBND_OP(log, "[-205,2)", "NaN");
		UBND_OP(log, "0", "-Inf");
		UBND_OP(log, "(2,4]", "(0.6931471805599453094,1.3862943611198906188]");
		UBND_OP(log, "(0,1)", "(-Inf,0)");
		UBND_OP(log, "-Inf", "NaN");
		UBND_OP(log, "10.123", "2.31481006261661465494044472");
		UBND_OP(log, "[-Inf,72381)", "NaN");
		UBND_OP(log, "(0,832.432]", "(-Inf,6.7243515368367869810983783837)");
		UBND_OP(log, "(4.62,21.32]", "(1.53039470509364748194379,3.059645599297643790718]");

		UBND_OP(exp, "[-20,2)", "[2.0611536224385578279659403e-9,7.3890560989306502272304)");
		UBND_OP(exp, "(2.5,4]", "(12.1824939607034734380701759,54.59815003314423907811]");
		UBND_OP(exp, "2.1", "8.166169912567650073449727410478631285183");
		UBND_OP(exp, "0.5509", "1.7348136477615857618355761641");
		UBND_OP(exp, "[-Inf,Inf]", "[0,Inf]");
		UBND_OP(exp, "(32.21,75]", "(9.7414871526611084e+13,3.733241996799001e+32]");
		UBND_OP(exp, "[-45.41,0.183]", "[1.8997111720014669e-20,1.200814408080830811743]");
		UBND_OP(exp, "[-92.32,-29)", "(8.0525500104695097629e-41,2.5436656473769229103033e-13]");

/* random log and exp test */
#define UBND_ROP(oper,op1) \
		ubnd_set_str(ub1, op1); \
		ubnd_##oper(ubr, ub1); \
		printf("%s ", #oper); ubnd_print(ub1); \
		printf(" = "); ubnd_print(ubr); putchar('\n')

		/* big number test */
		for (int i = 0; i < 0; i++){
			randomFloat(xstr,10,10);
			sign = rand() % 5;
			if (!sign) xstr[0] = '-';
			UBND_ROP(log,xstr);
			//UBND_ROP(exp,x);
		}

		/* random interval test */
		for (int i = 0; i < 5; i++){
			random_interval(xstr,max);
			UBND_ROP(log,xstr);
			//UBND_ROP(exp,interval);
		}

#define UBND_AOP(op1,oper,op2,chk) \
		ubnd_set_str(ub1, op1); \
		ubnd_set_str(ub2, op2); \
		ubnd_##oper(ubr, ub1, ub2); \
		ubnd_set_str(ubc, chk); \
		fail |= !(ok = ubnd_seq(ubr, ubc)); \
		printf("%s ", ok ? "OK  " : "FAIL"); \
		ubnd_print(ub1); printf(" %s ", #oper); ubnd_print(ub2); \
		printf(" = "); ubnd_print(ubr); putchar('\n')
		UBND_AOP("[2,6.5]", pow, "3.0", "[8,274.625]");
		UBND_AOP("3.2768e4", pow, "3", "35184372088832");
		UBND_AOP("[-1,2)", pow, "2", "[0,4)");
		UBND_AOP("(-2,2]", pow, "2", "[0,4]");
		UBND_AOP("(0.5,0.515625)", pow, "3", "(0.125,0.137088775634765625)");
		UBND_AOP("(-Inf,1.25]", pow, "2", "[0,Inf)");
		UBND_AOP("(-Inf,1.25]", pow, "3", "(-Inf,1.953125]");
		UBND_AOP("(-2,2]", pow, "3", "(-8,8]");
		UBND_AOP("2", pow, "0.5", "1.4142135623730950488");
		UBND_AOP("2.45", pow, "1.5", "(3.834808349609375,3.83489990234375)");
		UBND_AOP("(1,1.5625]", pow, "0.5", "(1,1.25]");
		UBND_AOP("2.44140625", pow, "0.25", "1.25");
		UBND_AOP("-2", pow, "2.1", "NaN");
		UBND_AOP("(-2,2)", pow, "2.1", "NaN");
		UBND_AOP("(-2,2)", pow, "(2,3)", "NaN");
		UBND_AOP("(0,2)", pow, "(2,3)", "(0,8)");
		UBND_AOP("(-2,2)", pow, "3", "(-8,8)");
		UBND_AOP("(-1,2)", pow, "2.0", "[0,4)");
		UBND_AOP("(-Inf,Inf)", pow, "2.0", "[0,Inf)");
		UBND_AOP("(-Inf,Inf)", pow, "(2.0,Inf)", "NaN");
		UBND_AOP("1.3", pow, "2.1", "1.73492633690415214192798878794939");

/* random power test */
#define UBND_RAOP(op1,oper,op2) \
		ubnd_set_str(ub1, op1); \
		ubnd_set_str(ub2, op2); \
		ubnd_##oper(ubr, ub1, ub2); \
		ubnd_print(ub1); printf(" %s ", #oper); ubnd_print(ub2); \
		printf(" = "); ubnd_print(ubr); putchar('\n')

		/* big number test */
		for (int i = 0; i < 0; i++){
			randomFloat(xstr,10,10);
			randomFloat(ystr,10,10);
			xsign = rand() % 5;
			if (!xsign) xstr[0] = '-';
			ysign = rand() % 5;
			if (!ysign) ystr[0] = '-';
			UBND_RAOP(xstr,pow,ystr);
		}

		/* random interval test */
		for (int i = 0; i < 5; i++){
			random_interval(xstr,max);
			random_interval(ystr,max);
			UBND_RAOP(xstr,pow,ystr);
		}

		/*
		UBND_AAOP("4",pow,"0.5");
		UBND_AAOP("16",pow,"0.25");
		UBND_AAOP("2",pow,"2");
		UBND_AAOP("-125",pow,"(1,2)");
		UBND_AAOP("[2,6.5]", pow, "3.0");
		UBND_AAOP("3.2768e4", pow, "3");
		UBND_AAOP("[-1,2)", pow, "2");
		UBND_AAOP("(-2,2]", pow, "2");
		UBND_AAOP("(0.5,0.515625)", pow, "3");
		UBND_AAOP("(-Inf,1.25]", pow, "2");
		UBND_AAOP("(-Inf,1.25]", pow, "3");
		UBND_AAOP("(-2,2]", pow, "3");
		UBND_AAOP("2", pow, "0.5");
		UBND_AAOP("2.45", pow, "1.5");
		
		UBND_AAOP("(1,2]", pow, "-1");
		UBND_AAOP("(1,1.5625]", pow, "-1");
		UBND_AAOP("[1,1.5625]", pow, "0.5");
		UBND_AAOP("2.44140625", pow, "0.25");
		
		UBND_AAOP("-2", pow, "2.1");
		UBND_AAOP("(-2,2)", pow, "2.1");
		UBND_AAOP("(-2,2)", pow, "(2,3)");
		UBND_AAOP("(-Inf,Inf)", pow, "(2.0,Inf)");
		UBND_AAOP("(0,2)", pow, "(2,3)");
		UBND_AAOP("(-2,2)", pow, "3");
		UBND_AAOP("(-1,2)", pow, "2.0");
		UBND_AAOP("(-Inf,Inf)", pow, "2.0");
		UBND_AAOP("1.3", pow, "2.1");
		
		UBND_AAOP("2.44140625", pow, "0.25");
		UBND_AAOP("2.562890625e-03",pow,"0.25");
		UBND_AAOP("0.25", pow, "0.5");
		UBND_AAOP("1.25", pow, "1.5");
		UBND_AAOP("12.25", pow, "-1.5");
		*/

		// -1000
		/*
		UBND_AOP("-212.5",pow,"-101");
		UBND_AOP("-11.5",pow,"-102");
		UBND_AOP("-1.15",pow,"-102");
		UBND_AOP("-0.845",pow,"-151");
		UBND_AOP("0.763",pow,"-251");
		UBND_AOP("2.43",pow,"-101");
		UBND_AOP("42.43",pow,"-101");
		UBND_AOP("226.37",pow,"-111");
		*/
		// -100
		/*
		UBND_AOP("-121.5",pow,"-20");
		UBND_AOP("-73.12",pow,"-11");
		UBND_AOP("-4.953",pow,"-36");
		UBND_AOP("-0.495",pow,"-78");
		UBND_AOP("0.952",pow,"-25");
		UBND_AOP("2.123",pow,"-32");
		UBND_AOP("26.943",pow,"-15");
		UBND_AOP("365.12",pow,"-15");
		*/

		// -10
		/*
		UBND_AOP("-196.2",pow,"-7");
		UBND_AOP("-34.12",pow,"-4");
		UBND_AOP("-3.246",pow,"-3");
		UBND_AOP("-0.823",pow,"-6");
		UBND_AOP("0.93",pow,"-2");
		UBND_AOP("5.742",pow,"-9");
		UBND_AOP("48.31",pow,"-4");
		UBND_AOP("348.9",pow,"-5");
		*/

		// -1
		/*
		UBND_AOP("0.93",pow,"-0.123");
		UBND_AOP("5.742",pow,"-0.837");
		UBND_AOP("48.31",pow,"-0.4");
		UBND_AOP("348.9",pow,"-0.23");
		*/
		// 1
		/*
		UBND_AOP("0.573",pow,"0.93");
		UBND_AOP("3.427",pow,"0.0463");
		UBND_AOP("94.75",pow,"0.5");
		UBND_AOP("734.28",pow,"0.437");
		*/

		// 10
		/*
		UBND_AOP("-116.2",pow,"17");
		UBND_AOP("-34.12",pow,"13");
		UBND_AOP("-3.246",pow,"15");
		UBND_AOP("-0.823",pow,"29");
		UBND_AOP("0.93",pow,"26");
		UBND_AOP("5.742",pow,"39");
		UBND_AOP("48.31",pow,"21.132");
		UBND_AOP("348.9",pow,"15");
		*/
		// 100
		/*
		UBND_AOP("-101.5",pow,"13");
		UBND_AOP("-23.12",pow,"14");
		UBND_AOP("-4.953",pow,"43");
		UBND_AOP("-0.495",pow,"92");
		UBND_AOP("0.952",pow,"91");
		UBND_AOP("2.123",pow,"73");
		UBND_AOP("16.943",pow,"13");
		UBND_AOP("105.32",pow,"19");
		*/
		// 1000
		/*
		UBND_AOP("-101.5",pow,"101");
		UBND_AOP("-13.12",pow,"103");
		UBND_AOP("-1.853",pow,"122");
		UBND_AOP("-0.895",pow,"232");
		UBND_AOP("0.952",pow,"116.2");
		UBND_AOP("2.123",pow,"113.13");
		UBND_AOP("26.943",pow,"102.21");
		UBND_AOP("365.12",pow,"124.62");
		*/
		//UBND_AOP("0.25",pow,"0.5");
		//UBND_AOP("2.562890625e-03",pow,"0.25");
		//UBND_AOP("720.023",pow,"2");

		tfail |= fail;
	}
#endif

/* Exp and Log function in ubnd time test */
#if 0
	unum_set_env(3,4);
	{
		UBND_VAR(xubnd[ITERS]);
		float xlist[ITERS];
		int ess, fss;
		int sign;
		float x;
		float max = 200000.0;

		unum_get_env(&ess, &fss);
		printf("\n# test time, env:%d,%d #\n", ess, fss);

		/* log function time test */
		#if 1
		{
			for (int i = 0; i <ITERS;i++){
				x = random_float(max);
				/* generate a random sign */
				sign = rand() % 5;
				if (!sign) x = -x; 
				ubnd_set_d(xubnd[i], x);
				xlist[i] = x;
				//printf("x: %f\n",x);
			}
		}
		#endif
		
		/* exp function time test */
		#if 0
		{
			for (int i = 0; i <ITERS;i++){
				x = random_float(max);
				sign = rand() % 5;
				if (!sign) x = -x;
				x = log(x); 
				ubnd_set_d(xubnd[i], x);
				xlist[i] = x;
				//printf("x: %f\n",x);
			}
		}
		#endif
		
		int64_t t1 = TimeInMicros();
		for (int iter = 0; iter < ITERS; iter++) {
			/* Unum time test */
			//ubnd_log(xubnd[iter], xubnd[iter]);
			//ubnd_exp(xubnd[iter], xubnd[iter]);

			/* IEEE time test */
			//xlist[iter] = log(xlist[iter]);
			//xlist[iter] = exp(xlist[iter]);
		}
		
		int64_t t2 = TimeInMicros();
    	printf("total time: %g ms; time per iteration: %g ms\n", (t2-t1)/1e3, (t2-t1)/1e3/ITERS);
	}
#endif

/* Pow function in ubnd time test */
#if 0
	unum_set_env(3,4);
	{
		int ess, fss;
		
		UBND_VAR(xubnd[ITERS]);
		UBND_VAR(yubnd[ITERS]);
		float x,y;
		float max = 30.0;
		int xsign, ysign;
		float xlist[ITERS], ylist[ITERS];
		
		unum_get_env(&ess, &fss);
		printf("\n# test time, env:%d,%d #\n", ess, fss);

		for (int i = 0; i <ITERS;i++){
			x = random_float(max);
			y = random_float(max);
			xsign = rand() % 3;
			ysign = rand() % 3;
			if (!xsign) x = -x; 
			if (!ysign) y = -y; 
			ubnd_set_d(xubnd[i], x);
			ubnd_set_d(yubnd[i], y);
			xlist[i] = x;
			ylist[i] = y;
			//printf("x: %f, y: %f\n", x, y);
		}
		
		
		int64_t t1 = TimeInMicros();
		for (int iter = 0; iter < ITERS; iter++) {
			//ubnd_pow(xubnd[ITERS], xubnd[iter],yubnd[iter]);
			//xlist[iter] = pow(xlist[iter],ylist[iter]);
		}
		int64_t t2 = TimeInMicros();
    	printf("total time: %g ms; time per iteration: %g ms\n", (t2-t1)/1e3, (t2-t1)/1e3/ITERS);
	}
#endif



	unum_clear_env();
	printf("\ntest %s\n", tfail ? "FAIL" : "OK");
	return(tfail ? EXIT_FAILURE : EXIT_SUCCESS);


}
