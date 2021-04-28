// [[Rcpp::plugins(openmp)]]
#define STRICT_R_HEADERS
#include <float.h>
#include <Rcpp.h>
#include <math.h>
#include <exception>

//OpenMP is not supported for macOS since R 4.0.0
#ifdef _OPENMP
#include <omp.h>
#endif // !_OPENMP

using namespace Rcpp;

struct Vgm
{
	int model;
	double nugget, psill, range;
};

int g_numOfIdsX, g_numOfIdsY;
DataFrame g_discretePointsX, g_discretePointsY;
Vgm g_ptVgmModelX, g_ptVgmModelY;
NumericMatrix g_areaDistByCentroidX, g_areaDistByCentroidXY;
std::vector<NumericMatrix> g_areaDistByPtsX, g_areaDistByPtsY, g_areaDistByPtsXY;
std::vector<NumericMatrix> g_areaWeightByPtsX, g_areaWeightByPtsY, g_areaWeightByPtsXY;


inline NumericVector CalcVariogramSimple(const Vgm& vgm, const NumericVector& dist, bool bCov = false)
{
	int nSize = dist.size();
	NumericVector semivar(nSize);
	switch (vgm.model)
	{
	case 1: // Exp
		// semivar = vgm.nugget + vgm.psill * (1 - exp(-dist / vgm.range));
			//omp_set_num_threads(2);
		for (int i = 0; i < nSize; i++)
			semivar[i] = vgm.nugget + vgm.psill * (1 - exp(-dist[i] / vgm.range));
		break;
	case 2: // Gau
		//semivar = vgm.nugget + vgm.psill * (1 - exp(-pow(dist, 2) / (vgm.range * vgm.range)));
		for (int i = 0; i < nSize; i++)
			semivar = vgm.nugget + vgm.psill * (1 - exp(-(dist[i] * dist[i]) / (vgm.range * vgm.range)));
		break;
	case 3: //Sph
		//semivar = vgm.nugget + vgm.psill * (1.5 * dist / vgm.range - 0.5 * pow(dist / vgm.range, 3));
		//semivar[dist >= vgm.range] = vgm.nugget + vgm.psill;
		for (int i = 0; i < nSize; i++)
		{
			semivar[i] = vgm.nugget + vgm.psill * (1.5 * dist[i] / vgm.range - 0.5 * pow(dist[i] / vgm.range, 3));
			if (dist[i] >= vgm.range) semivar[i] = vgm.nugget + vgm.psill;
		}
		break;
	default:
		return NumericVector::create();
		break;
	}

	if (bCov) semivar = (vgm.nugget + vgm.psill) - semivar;
	if (Rf_isMatrix(dist)) semivar.attr("dim") = dist.attr("dim");
	return semivar;
}


inline NumericVector CalcVariogram_gstat(const DataFrame& vgm, const NumericVector& dist, bool bCov = false)
{
	Rcpp::Environment gstat("package:gstat");
	Rcpp::Function variogramLine = gstat["variogramLine"];
	NumericVector mSvar = variogramLine(vgm, Named("dist_vector") = dist, Named("covariance") = bCov);
	return mSvar;
}


inline double CalcWeightedVariogram(const Vgm& vgm, const NumericVector& dist, const NumericVector& weight)
{
	int nSize = dist.size();
	double semivar = 0.0;

	switch (vgm.model)
	{
	case 1: // Exp
		for (int i = 0; i < nSize; i++)
			semivar += weight[i] * (vgm.nugget + vgm.psill * (1 - exp(-dist[i] / vgm.range)));
		break;
	case 2: // Gau
		for (int i = 0; i < nSize; i++)
			semivar += weight[i] * (vgm.nugget + vgm.psill * (1 - exp(-(dist[i] * dist[i]) / (vgm.range * vgm.range))));
		break;
	case 3: //Sph
		double temp;
		for (int i = 0; i < nSize; i++)
		{
			temp = vgm.nugget + vgm.psill * (1.5 * dist[i] / vgm.range - 0.5 * pow(dist[i] / vgm.range, 3));
			if (dist[i] >= vgm.range) temp = vgm.nugget + vgm.psill;
			semivar += weight[i] * temp;
		}
		break;
	default:
		return NA_REAL;
		break;
	}

	return semivar;
}


Vgm VgmFromDf(const DataFrame& vgmDf) {
	Vgm vgm;
	String modelName;
	if (vgmDf.nrows() == 1) { // no nugget
		modelName = as<CharacterVector>(vgmDf[0])[0];
		vgm.nugget = 0.0;
		vgm.psill = as<NumericVector>(vgmDf[1])[0];
		vgm.range = as<NumericVector>(vgmDf[2])[0];
	}
	else {
		modelName = as<CharacterVector>(vgmDf[0])[1];
		vgm.nugget = as<NumericVector>(vgmDf[1])[0];
		vgm.psill = as<NumericVector>(vgmDf[1])[1];
		vgm.range = as<NumericVector>(vgmDf[2])[1];
	}

	if (modelName == "Exp")
		vgm.model = 1;
	else if (modelName == "Gau")
		vgm.model = 2;
	else if (modelName == "Sph")
		vgm.model = 3;
	else
		vgm.model = -1;

	return vgm;
}


// [[Rcpp::export]]
RObject variogramLineSimple(const DataFrame& vgmModel, const NumericVector& dist, bool bCov = false)
{
	Vgm vgm = VgmFromDf(vgmModel);
	NumericVector semivar = CalcVariogramSimple(vgm, dist, bCov);
	if (bCov || Rf_isMatrix(dist))
		return semivar;
	else
		return DataFrame::create(Named("dist") = dist, Named("gamma") = semivar);
}


inline LogicalVector CompareCharacter(CharacterVector& vec, String scalar)
{
	LogicalVector indx(vec.size());
	for (int i = 0; i < vec.size(); i++)
		indx[i] = vec[i] == scalar;
	return indx;
}


/*  The function is modified from sp package (Copyright by Roger Bivand (C) 2005-2009), thanks Roger Bivand.  */
# define POWDI(x,i) pow(x,i)
double sp_gcdist(double lon1, double lon2, double lat1, double lat2) {

	double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C;
	double w, R, a, f, D, H1, H2;
	double lat1R, lat2R, lon1R, lon2R, DE2RA;
	double dist;

	DE2RA = M_PI / 180;
	a = 6378.137;              /* WGS-84 equatorial radius in km */
	f = 1.0 / 298.257223563;     /* WGS-84 ellipsoid flattening factor */

	if (fabs(lat1 - lat2) < DBL_EPSILON) {
		if (fabs(lon1 - lon2) < DBL_EPSILON) {
			dist = 0.0;
			return dist;
			/* Wouter Buytaert bug caught 100211 */
		}
		else if (fabs((fabs(lon1) + fabs(lon2)) - 360.0) < DBL_EPSILON) {
			dist = 0.0;
			return dist;
		}
	}
	lat1R = lat1 * DE2RA;
	lat2R = lat2 * DE2RA;
	lon1R = lon1 * DE2RA;
	lon2R = lon2 * DE2RA;

	F = (lat1R + lat2R) / 2.0;
	G = (lat1R - lat2R) / 2.0;
	L = (lon1R - lon2R) / 2.0;

	/*
	printf("%g %g %g %g; %g %g %g\n",  *lon1, *lon2, *lat1, *lat2, F, G, L);
	*/

	sinG2 = POWDI(sin(G), 2);
	cosG2 = POWDI(cos(G), 2);
	sinF2 = POWDI(sin(F), 2);
	cosF2 = POWDI(cos(F), 2);
	sinL2 = POWDI(sin(L), 2);
	cosL2 = POWDI(cos(L), 2);

	S = sinG2 * cosL2 + cosF2 * sinL2;
	C = cosG2 * cosL2 + sinF2 * sinL2;

	w = atan(sqrt(S / C));
	R = sqrt(S * C) / w;

	D = 2 * w * a;
	H1 = (3 * R - 1) / (2 * C);
	H2 = (3 * R + 1) / (2 * S);

	dist = D * (1 + f * H1 * sinF2 * cosG2 - f * H2 * cosF2 * sinG2);
	return dist;
}


// [[Rcpp::export]]
NumericMatrix spDistsNN(NumericVector& x1, NumericVector& y1, NumericVector& x2, NumericVector& y2, bool longlat = false)
{
	int n = x1.size();
	int m = x2.size();
	NumericMatrix d(n, m);

	if (!longlat) {
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				d(i, j) = hypot(x1[i] - x2[j], y1[i] - y2[j]);
	}
	else {
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				d(i, j) = sp_gcdist(x1[i], x2[j], y1[i], y2[j]);
	}

	return d;
}


// [[Rcpp::export]]
NumericMatrix outerProd(const NumericVector& v1, const NumericVector& v2) {
	NumericMatrix r(v1.size(), v2.size());
	for (int i = 0; i < v1.size(); i++) {
		for (int j = 0; j < v2.size(); j++) {
			r(i, j) = v1[i] * v2[j];
		}
	}
	return r;
}


// [[Rcpp::export]]
void svAreaCloudByPointVgmInit(const DataFrame& discretePoints, const NumericMatrix& areaDistByCentroid, bool longlat = false)
{
	g_discretePointsX = discretePoints;
	g_areaDistByCentroidX = areaDistByCentroid;

	CharacterVector areaId = as<CharacterVector>(g_discretePointsX[0]);
	g_numOfIdsX = sort_unique(areaId).size();

	//distances between discretized area points
	g_areaDistByPtsX.clear();
	g_areaWeightByPtsX.clear();
	/*
	for (int i = 0; i < g_numOfIdsX - 1; i++) {
		g_areaDistByPtsX.push_back(as<NumericMatrix>(as<List>(areaDistByPts[i])[i]));
		g_areaWeightByPtsX.push_back(as<NumericMatrix>(as<List>(areaWeightByPts[i])[i]));

		for (int j = i + 1; j < g_numOfIdsX; j++) {
			g_areaDistByPtsX.push_back(as<NumericMatrix>(as<List>(areaDistByPts[j])[j]));
			g_areaWeightByPtsX.push_back(as<NumericMatrix>(as<List>(areaWeightByPts[j])[j]));

			g_areaDistByPtsX.push_back(as<NumericMatrix>(as<List>(areaDistByPts[i])[j]));
			g_areaWeightByPtsX.push_back(as<NumericMatrix>(as<List>(areaWeightByPts[i])[j]));
		}
	}
	*/

	CharacterVector uId = sort_unique(areaId);
	NumericVector x1, y1, x2, y2, w1, w2;
	LogicalVector areaIndexI, areaIndexJ;

	//int N = g_numOfIdsX;
	//int indexI = 0, indexJ = 0;
	int i, j;

	/*
	g_areaDistByPtsX.reserve((N - 1) * (N + 1));
	g_areaWeightByPtsX.reserve((N - 1) * (N + 1));
	for (int i = 0; i < (N-1) * (N+1); i++)
	{
		g_areaDistByPtsX.push_back(NumericMatrix(2, 2));
		g_areaWeightByPtsX.push_back(NumericMatrix(2, 2));
	}
	*/

//#pragma omp parallel for private(i, j, areaIndexI, areaIndexJ, indexI, indexJ, x1, y1, x2, y2, w1, w2)
	for (i = 0; i < uId.size() - 1; i++)
	{
		areaIndexI = CompareCharacter(areaId, uId[i]);
		x1 = as<NumericVector>(discretePoints[1])[areaIndexI];
		y1 = as<NumericVector>(discretePoints[2])[areaIndexI];
		w1 = as<NumericVector>(discretePoints[3])[areaIndexI];

		g_areaDistByPtsX.push_back(spDistsNN(x1, y1, x1, y1));
		g_areaWeightByPtsX.push_back(outerProd(w1, w1));

		/*
		indexI = (N * (N + 1) - (N - i) * (N + 1 - i)) - i;
		g_areaDistByPtsX[indexI] = spDists(x1, y1, x1, y1);
		g_areaWeightByPtsX[indexI] = outer(w1, w1);
		*/
		for (j = i+1; j < uId.size(); j++)
		{
			areaIndexJ = CompareCharacter(areaId, uId[j]);
			x2 = as<NumericVector>(discretePoints[1])[areaIndexJ];
			y2 = as<NumericVector>(discretePoints[2])[areaIndexJ];
			w2 = as<NumericVector>(discretePoints[3])[areaIndexJ];

			g_areaDistByPtsX.push_back(spDistsNN(x2, y2, x2, y2));
			g_areaDistByPtsX.push_back(spDistsNN(x1, y1, x2, y2));
			g_areaWeightByPtsX.push_back(outerProd(w2, w2));
			g_areaWeightByPtsX.push_back(outerProd(w1, w2));

			/*
			indexJ = indexI + (j - i) * 2 - 1;
			g_areaDistByPtsX[indexJ] = spDists(x2, y2, x2, y2);
			g_areaDistByPtsX[indexJ + 1] = spDists(x1, y1, x2, y2);
			g_areaWeightByPtsX[indexJ] = outer(w2, w2);
			g_areaWeightByPtsX[indexJ + 1] = outer(w1, w2);
			*/
		}
	}
}


// [[Rcpp::export]]
void svAreaCloudByPointVgmEnd()
{
	g_areaDistByPtsX.clear();
	g_areaWeightByPtsX.clear();
}


// test
DataFrame svAreaCloudByPointVgm_gstat(const DataFrame& ptVgmModel)
{
	Rcpp::Environment gstat("package:gstat");
	Rcpp::Function variogramLine = gstat["variogramLine"];

	NumericMatrix dg(g_numOfIdsX * (g_numOfIdsX - 1) / 2, 2);
	NumericVector mSvar;
	double g11, g12, g22, g;

	int N = g_numOfIdsX;
	int indexI = 0, indexJ = 0, indexK = 0;
	int i, j;

//#pragma omp parallel for private(i, j, indexI, indexJ, indexK, g11, g12, g22, g)
	for (i = 0; i <= g_numOfIdsX - 2; i++) {
		indexI = (N * (N + 1) - (N - i) * (N + 1 - i)) - i;
		mSvar = variogramLine(ptVgmModel, Named("dist_vector") = g_areaDistByPtsX[indexI]);
		g11 = sum(g_areaWeightByPtsX[indexI] * mSvar);

		for (j = i + 1; j <= g_numOfIdsX - 1; j++) {
			indexJ = indexI + (j - i) * 2 - 1;
			mSvar = variogramLine(ptVgmModel, Named("dist_vector") = g_areaDistByPtsX[indexJ]);
			g22 = sum(g_areaWeightByPtsX[indexJ] * mSvar);

			mSvar = variogramLine(ptVgmModel, Named("dist_vector") = g_areaDistByPtsX[indexJ + 1]);
			g12 = sum(g_areaWeightByPtsX[indexJ + 1] * mSvar);

			g = g12 - (g11 + g22) / 2.0;
			indexK = ((N * (N + 1) - (N - i) * (N + 1 - i))) / 2 + (j - i) - (i + 1);
			dg(indexK, 0) = g_areaDistByCentroidX(i, j);
			dg(indexK, 1) = g;
		}
	}

	DataFrame df = as<DataFrame>(dg);
	df.names() = CharacterVector::create("dist", "gamma");
	return df;
}


// [[Rcpp::export]]
void ataSetNumberOfThreadsForOMP(int num)
{
#ifdef _OPENMP
	int n = omp_get_num_procs();
	if(num > 0 && num <= n)
		omp_set_num_threads(num);
#else
  Rcout << "Not supported since OPENMP is not available!\n";
#endif
}


// [[Rcpp::export]]
DataFrame svAreaCloudByPointVgm(const DataFrame& ptVgmModel)
{
	// Rprintf("here!");
	NumericMatrix dg(g_numOfIdsX * (g_numOfIdsX - 1) / 2, 2);
	Vgm vgm = VgmFromDf(ptVgmModel);

	NumericVector mSvar;
	double g11, g12, g22, g;

	int N = g_numOfIdsX;
	int indexI = 0, indexJ = 0, indexK = 0;
	int i, j;

#pragma omp parallel for private(i, j, indexI, indexJ, indexK, g11, g12, g22, g)
	for (i = 0; i <= g_numOfIdsX - 2; i++) {
		indexI = (N * (N + 1) - (N - i) * (N + 1 - i)) - i;
		g11 = CalcWeightedVariogram(vgm, g_areaDistByPtsX[indexI], g_areaWeightByPtsX[indexI]);

		for (j = i + 1; j <= g_numOfIdsX - 1; j++) {
			indexJ = indexI + (j - i) * 2 - 1;
			g22 = CalcWeightedVariogram(vgm, g_areaDistByPtsX[indexJ], g_areaWeightByPtsX[indexJ]);
			g12 = CalcWeightedVariogram(vgm, g_areaDistByPtsX[indexJ + 1], g_areaWeightByPtsX[indexJ + 1]);
			g = g12 - (g11 + g22) / 2.0;

			indexK = ((N * (N + 1) - (N - i) * (N + 1 - i))) / 2 + (j - i) - (i + 1);
			dg(indexK, 0) = g_areaDistByCentroidX(i, j);
			dg(indexK, 1) = g;
		}
	}

	DataFrame df = as<DataFrame>(dg);
	df.names() = CharacterVector::create("dist", "gamma");
	return df;
}


// [[Rcpp::export]]
void crossSvAreaCloudByPointVgmInit(
	const DataFrame& discretePointsX, const DataFrame& discretePointsY,
	const DataFrame& ptVgmModelX, const DataFrame& ptVgmModelY,
	const NumericMatrix& areaDistByCentroidXY,
	const List& areaDistByPtsX, const List& areaDistByPtsY, const List& areaDistByPtsXY,
	const List& areaWeightByPtsX, const List& areaWeightByPtsY, const List& areaWeightByPtsXY)
{
	g_discretePointsX = discretePointsX;
	g_discretePointsY = discretePointsY;
	g_areaDistByCentroidXY = areaDistByCentroidXY;

	g_ptVgmModelX = VgmFromDf(ptVgmModelX);
	g_ptVgmModelY = VgmFromDf(ptVgmModelY);

	g_numOfIdsX = sort_unique(as<CharacterVector>(g_discretePointsX[0])).size();
	g_numOfIdsY = sort_unique(as<CharacterVector>(g_discretePointsY[0])).size();

	g_areaDistByPtsX.clear();
	g_areaDistByPtsY.clear();
	g_areaDistByPtsXY.clear();
	g_areaWeightByPtsX.clear();
	g_areaWeightByPtsY.clear();
	g_areaWeightByPtsXY.clear();
	for (int i = 0; i < g_numOfIdsX; i++)
	{
		g_areaDistByPtsX.push_back(as<NumericMatrix>(areaDistByPtsX[i]));
		g_areaWeightByPtsX.push_back(as<NumericMatrix>(areaWeightByPtsX[i]));
		for (int j = 0; j < g_numOfIdsY; j++)
		{
			g_areaDistByPtsXY.push_back(as<NumericMatrix>(as<List>(areaDistByPtsXY[i])[j]));
			g_areaWeightByPtsXY.push_back(as<NumericMatrix>(as<List>(areaWeightByPtsXY[i])[j]));
		}
	}
	for (int j = 0; j < g_numOfIdsY; j++)
	{
		g_areaDistByPtsY.push_back(as<NumericMatrix>(areaDistByPtsY[j]));
		g_areaWeightByPtsY.push_back(as<NumericMatrix>(areaWeightByPtsY[j]));
	}
}


// [[Rcpp::export]]
DataFrame crossSvAreaCloudByPointVgm(const DataFrame& xyPointCrossVgm)
{
	NumericMatrix dg(g_numOfIdsX * g_numOfIdsX, 2);
	Vgm vgmXY = VgmFromDf(xyPointCrossVgm);

	NumericVector mSvar;
	double g11, g12, g22, g;
	int i, j, indexJ;

#pragma omp parallel for private(i, j, indexJ, g11, g12, g22, g)
	for (i = 0; i < g_numOfIdsX; i++) {
		g11 = CalcWeightedVariogram(g_ptVgmModelX, g_areaDistByPtsX[i], g_areaWeightByPtsX[i]);

		for (j = 0; j < g_numOfIdsY; j++) {
			g22 = CalcWeightedVariogram(g_ptVgmModelY, g_areaDistByPtsY[j], g_areaWeightByPtsY[j]);

			indexJ = i * g_numOfIdsY + j;
			g12 = CalcWeightedVariogram(vgmXY, g_areaDistByPtsXY[indexJ], g_areaWeightByPtsXY[indexJ]);
			g = g12 - (g11 + g22) / 2.0;

			dg(indexJ, 0) = g_areaDistByCentroidXY(i, j);
			dg(indexJ, 1) = g;
		}
	}

	DataFrame df = as<DataFrame>(dg);
	df.names() = CharacterVector::create("dist", "gamma");
	return df;
}


// [[Rcpp::export]]
void crossSvAreaCloudByPointVgmEnd()
{
	g_areaDistByPtsX.clear();
	g_areaDistByPtsY.clear();
	g_areaDistByPtsXY.clear();
	g_areaWeightByPtsX.clear();
	g_areaWeightByPtsY.clear();
	g_areaWeightByPtsXY.clear();
}

