/*
 * transformation.cpp
 *
 *  Created on: Apr 24, 2012
 *      Author: wjiang2
 */

#include "include/transformation.hpp"

PARAM_VEC::iterator findTransFlag(PARAM_VEC & pVec, string name){
	PARAM_VEC::iterator it;
	for(it=pVec.begin();it!=pVec.end();it++)
	{
		if(it->param.compare(name)==0)
			break;
	}
	return it;
}
trans_map trans_local::cloneTransMap(){

	trans_map res;
	/*
	 * clone trans map
	 */

	for(trans_map::iterator it=transformations.begin();it!=transformations.end();it++)
	{
		transformation * curTran=it->second;
		if(curTran!=NULL)
		{
			cout<<"cloning transformatioin:"<<curTran->channel<<endl;
			res[it->first]=curTran->clone();
		}
	}
	return res;
}


/*
 * transformation
 */



transformation::transformation(){
	/*
	 * if it is pure transformation object,then assume calibration is directly read from ws
	 * so there is no need to compute calibration
	 */
	isComputed=true;
	isGateOnly=false;

}

double mylog10(double x) {
	return x>0?log10(x):0;
}
/*
 * these transforming functions change the input data
 */

void logTrans::transforming(valarray<double> & input){

		/*
		 * clean the data before log
		 */
		input=input.apply(mylog10);
//		input=log10(input);

}
void linTrans::transforming(valarray<double> & input){

		input*=64;
}
void transformation::transforming(valarray<double> & input){
		if(!calTbl.isInterpolated)
			throw(domain_error("calibration table not interpolated yet!"));
		input=calTbl.transforming(input);
}

transformation * trans_local::getTran(string channel){
	transformation * res;
	if(channel.compare("Time")==0||channel.compare("time")==0)
		res=NULL;


	trans_map::iterator it=transformations.find(channel);
	if(it==transformations.end())
		res=NULL;
	else
		res=it->second;

	return res;
}
/*
 *biexpTrans
 */
biexpTrans * biexpTrans::clone(){

	return new biexpTrans(*this);
}


biexpTrans::biexpTrans(){
	channelRange=4096;
	maxValue=262144;
	pos=4.5;
	neg=0;
	widthBasis=-10;
	isComputed=false;

}
/*
 * directly translated from java routine from tree star
 */
 double logRoot(double b, double w)
{
	double xLo = 0;
	double xHi = b;
	double d = (xLo + xHi) / 2;
	double dX = abs((long) (xLo - xHi));
	double dXLast = dX;
	double fB = -2 * log(b) + w * b;
	double f = 2. * log(d) + w * b + fB;
	double dF = 2 / d + w;
	if (w == 0) return b;
	for (long i = 0; i < 100; i++)
	{
		if (((d - xHi) * dF - f) * ((d - xLo) * dF - f) >= 0 ||
				abs((long) (2 * f)) > abs((long) (dXLast * dF)))
		{
			dX = (xHi - xLo) / 2;
			d = xLo + dX;
			if (d == xLo)
				return d;
		}
		else
		{
			dX = f / dF;
			double t = d;
			d -= dX;
			if (d == t)
				return d;
		}
		if (abs((long) dX) < 1.0e-12)
			return d;
		dXLast = dX;
		f = 2 * log(d) + w * d + fB;
		dF = 2 / d + w;
		if (f < 0) xLo = d;
		else xHi = d;
	}
	return d;
}
 /*
  * directly translated from java routine from tree star
  */

void logInterpolate(double *f, int i, long n, double x)
{
	double minVal = f[i];
	double maxVal = x;
	double logMinVal = log(minVal);
	double logMaxVal = log(maxVal);
	for (int j = i; j < n; j++)
	{
		float frxn = (float)(j - i) / (float)(n - i);
		double curVal = frxn * (logMaxVal - logMinVal) + logMinVal;
		f[j] = exp(curVal);
	}
}
void biexpTrans::computCalTbl(){
	/*
	 * directly translated from java routine from tree star
	 */

	double ln10 = log(10.0);
	double decades = pos;
	double lowScale = widthBasis;
	double width = log10(-lowScale);

	if (width < 0.5 || width > 3) width = 0.5;
	decades -= width / 2;
	double extra = neg;
	if (extra < 0) extra = 0;
	extra += width / 2;

	int zeroChan = (int)(extra * channelRange / (extra + decades));
	zeroChan = min(zeroChan, channelRange / 2);

	if (zeroChan > 0) decades = extra * channelRange / zeroChan;
	width /= 2 * decades;        // 1.1

	double maximum = maxValue;
	double positiveRange = ln10 * decades;
	double minimum = maximum / exp(positiveRange);
	double negativeRange = logRoot(positiveRange, width);

	double *positive = new double[channelRange + 1];
	double *negative = new double[channelRange + 1];
	positive[0] = negative[0] = 1.;
	logInterpolate(positive, 0, channelRange + 1, exp(positiveRange));
	logInterpolate(negative, 0, channelRange + 1, exp(-negativeRange));

	double s = exp((positiveRange + negativeRange) * (width + extra / decades));
	int j;
	for (j = 0; j < channelRange + 1; j++)
		negative[j] *= s;
	s = positive[zeroChan] - negative[zeroChan];
	for (j = zeroChan; j < channelRange + 1; j++)
		positive[j] = minimum * (positive[j] - negative[j] - s);
	for (j = 0; j < zeroChan; j++)
		positive[j] = -positive[2 * zeroChan - j];

	/*
	 * save the calibration table
	 */
	calTbl.caltype="flowJo";
	calTbl.spline_method=2;
	calTbl.init(channelRange+1);

	for (int chan = 0; chan <= channelRange; chan++)
	{
		calTbl.y[chan] =chan;
		calTbl.x[chan] = positive[chan];
	}

	isComputed=true;

	delete[] positive;
	delete[] negative;

}
/*
 * calTrans
 */
//calTrans * calTrans::clone(){
//
//	return new calTrans(*this);
//}
