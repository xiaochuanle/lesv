#include "lof_point.hpp"
#include "lof.hpp"

#include <iomanip>
#include <iostream>
#include <initializer_list>

using namespace std;
using namespace hbn_lof;

#if 1
typedef LofPoint<2> inst_type;

static initializer_list<inst_type> test_instances1 = {
 inst_type(-4.8447532242074978, -5.6869538132901658),
 inst_type(1.7265577109364076, -2.5446963280374302),
 inst_type(-1.9885982441038819, 1.705719643962865),
 inst_type(-1.999050026772494, -4.0367551415711844),
 inst_type(-2.0550860126898964, -3.6247409893236426),
 inst_type(-1.4456945632547327, -3.7669258809535102),
 inst_type(-4.6676062022635554, 1.4925324371089148),
 inst_type(-3.6526420667796877, -3.5582661345085662),
 inst_type(6.4551493172954029, -0.45434966683144573),
 inst_type(-0.56730591589443669, -5.5859532963153349),
 inst_type(-5.1400897823762239, -1.3359248994019064),
 inst_type(5.2586932439960243, 0.032431285797532586),
 inst_type(6.3610915734502838, -0.99059648246991894),
 inst_type(-0.31086913190231447, -2.8352818694180644),
 inst_type(1.2288582719783967, -1.1362795178325829),
 inst_type(-0.17986204466346614, -0.32813130288006365),
 inst_type(2.2532002509929216, -0.5142311840491649),
 inst_type(-0.75397166138399296, 2.2465141276038754),
 inst_type(1.9382517648161239, -1.7276112460593251),
 inst_type(1.6809250808549676, -2.3433636210337503),
 inst_type(0.68466572523884783, 1.4374914487477481),
 inst_type(2.0032364431791514, -2.9191062023123635),
 inst_type(-1.7565895138024741, 0.96995712544043267),
 inst_type(3.3809644295064505, 6.7497121359292684),
 inst_type(-4.2764152718650896, 5.6551328734397766),
 inst_type(-3.6347215445083019, -0.85149861984875741),
 inst_type(-5.6249411288060385, -3.9251965527768755),
 inst_type(4.6033708001912093, 1.3375110154658127),
 inst_type(-0.685421751407983, -0.73115552984211407),
 inst_type(-2.3744241805625044, 1.3443896265777866)
};


int main()
{
    vector<inst_type> instances(test_instances1);
    vector<LofOutlier<2>> outliers = GetOutliers(5, instances);
    cout << "number of outlier: " << outliers.size() << endl;
	for (int i = 0; i < outliers.size(); i++)
	{
		cout << setiosflags(ios::fixed) << setprecision(20) << outliers.at(i).lof << " ( ";
		cout << setiosflags(ios::fixed) << setprecision(20) << outliers.at(i).instance[0]<< " , ";
		cout << setiosflags(ios::fixed) << setprecision(20) << outliers.at(i).instance[1]<< " ) "<<endl;
	}
}
#endif