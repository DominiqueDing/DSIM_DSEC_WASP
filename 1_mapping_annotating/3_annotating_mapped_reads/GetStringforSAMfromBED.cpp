//Usage: ./GetStringforSAMfomBED <BED> <incomplete.SAM> <complete.SAM>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream> 
#include <algorithm> //std::sort
#include <functional> //std::greater

using namespace std;

int main(int argc, char const *argv[])
{
	string line;
	//import BED:
	vector<string> BED_v;
	ifstream BED (argv[1]);
	while (getline(BED,line))
	{
		if (BED.good())
		{
			int count = 0; size_t found = 0, found1, found2, found3, found4;
			while (count != 14)
			{
				found = line.find('\t', found+1);
				count ++;
				if (count == 3)
				{
					found1 = found;
				}
				else if (count == 4)
				{
					found2 = found;
				}
				else if (count == 11)
				{
					found3 = found;
				}
				else if (count == 13)
				{
					found4 = found;
				}
			}
			string line_tmp = line.substr(found1+1,found2-found1) + line.substr(found3,found4-found3);
			line.assign(line_tmp);//extracting only useful bits: seqID seqString qualString
			BED_v.push_back(line);
			line.clear();
		}
	}
	if (BED_v.empty())
	{
		cerr << "Something went wrong reading BED file!\n";
		return -1;
	}
	else cout << "Read in BED done! \n";

	sort(BED_v.begin(),BED_v.end());//sort BED_v
	if (BED_v.empty())
	{
		cerr << "Something went wrong sorting BED file!\n";
		return -1;
	}
	else if (adjacent_find(BED_v.begin(), BED_v.end(), greater<string>())!=BED_v.end())
	{
		cerr << "BED not sorted!\n";
		return -1;
	}
	else cout << "BED sorted!\n";
	
	
	//import incomplete SAM
	vector<string> in_SAM_v;
	ifstream inSAM (argv[2]); line.clear();
	while (getline(inSAM,line))
	{
		if (inSAM.good())
		{
			int count = 0; size_t found=0;
			while(count != 9)
			{
				found = line.find('\t',found+1);
				count ++;
			}
			line.erase(found,string::npos);//get rid of the last 2 fields
			in_SAM_v.push_back(line);
			line.clear();
		}
	}
	if (in_SAM_v.empty())
	{
		cerr << "Something went wrong reading SAM file!\n";
		return -1;
	}
	else cout << "Read in SAM done! \n";

	sort(in_SAM_v.begin(),in_SAM_v.end());//sort in_SAM_v
	if (in_SAM_v.empty())
	{
		cerr << "Something went wrong sorting SAM file!\n";
		return -1;
	}
	else if (adjacent_find(in_SAM_v.begin(), in_SAM_v.end(), greater<string>())!=in_SAM_v.end())
	{
		cerr << "SAM not sorted!\n";
		return -1;
	}
	else cout << "SAM sorted!\n";

	in_SAM_v.size() == BED_v.size();
	ofstream out_SAM (argv[3]);
	for (int i = 0; i < in_SAM_v.size(); ++i)
	{
		string samID = in_SAM_v[i].substr(0,in_SAM_v[i].find('\t'));
		string bedID = BED_v[i].substr(0,BED_v[i].find('\t'));
		if (samID == bedID)
		{
			BED_v[i].erase(0,in_SAM_v[i].find('\t')+1);
			out_SAM << in_SAM_v[i] << BED_v[i] << '\n';
		}
		else
		{
			cerr << "SAM and BED don't match!\n";
			return -1;
		}
	}
	cout << argv[3] << " complete!\n";


	return 0;
}