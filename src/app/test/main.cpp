#include "../../algo/hbn_traceback.h"
#include "../../algo/hbn_traceback_aux.h"
#include "../../algo/refine_align.h"
#include "../../algo/approx_semi_gapped_align.h"
#include "../../corelib/line_reader.h"
#include "../../corelib/khash.h"
#include "../../corelib/m4_record.h"
#include "../../corelib/fasta.h"
#include "../../corelib/seq_tag_report.h"
#include "../../corelib/string2hsp.h"
#include "../../corelib/ksort.h"
#include "../../corelib/cstr_util.h"
#include "../../corelib/seqdb.h"
#include "../../ncbi_blast/setup/hsp2string.h"
#include "../../corelib/cstr_util.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <ctype.h>
#include <pthread.h>

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

#if 0
int main(int argc, char* argv[])
{
	ifstream in(argv[1]);
	string line;
	vector<ncbi::CTempString> components;
	const ncbi::CTempString kDelim("\t");

	while(getline(in, line)) {
		components.clear();
		NStr::Split(line, kDelim, components);
		hbn_assert(components.size() == 12, "%zu", components.size());
		int ssize = NStr::StringToInt(components[11]);
		if (ssize != 7088241) continue;
		double perc_identity = NStr::StringToDouble(components[2]);
		//if (perc_identity < 94.0) continue;

		int qoff = NStr::StringToInt(components[5]);
		int qend = NStr::StringToInt(components[6]);
		int qsize = NStr::StringToInt(components[7]);
		const int E = 200;
		if (!(qoff <= E && qsize - qend <= E)) continue;

		int soff = NStr::StringToInt(components[9]);
		int send = NStr::StringToInt(components[10]);
		//cerr << soff << '\t' << send << '\t' << ssize << '\t' << perc_identity << endl;
		const int bp = 7088241;
		if (soff <= bp && send >= bp) cout << line << endl;
	}

	in.close();
	return 0;
}
#endif

#if 0
static void
cov_stat_for_one_contig(const char* m4_path, CSeqDB* sdb, int sid)
{
	int ssize = seqdb_seq_size(sdb, sid);
	if (ssize != 18072166) return;
	const char* sname = seqdb_seq_name(sdb, sid);
	HBN_LOG("cov stats for subject %s, length = %d", sname, ssize);
	int* cov_stats = (int*)calloc(ssize, sizeof(int));

	ifstream in(m4_path);
	string line;
	vector<ncbi::CTempString> components;
	const ncbi::CTempString kDelim("\t");

	while(getline(in, line)) {
		components.clear();
		NStr::Split(line, kDelim, components);
		hbn_assert(components.size() == 12, "%zu", components.size());

		string name = components[1];
		if (strcmp(sname, name.c_str())) continue;

		int qoff = NStr::StringToInt(components[5]);
		int qend = NStr::StringToInt(components[6]);
		int qsize = NStr::StringToInt(components[7]);
		int soff = NStr::StringToInt(components[9]);
		int send = NStr::StringToInt(components[10]);
		double perc_identity = NStr::StringToDouble(components[2]);

		const int E = 200;
		if (!(qoff <= E && qsize - qend <= E)) continue;
		if (perc_identity < 92.0) continue;

		if (soff < 10494019 && send > 10494019) cout << line << endl;

		for (int i = soff; i < send; ++i) ++cov_stats[i];
	}
	in.close();

	int i = 0;
	int cov = 5;
	while (i < ssize) {
		while (i < ssize && cov_stats[i] < cov) ++i;
		if (i >= ssize) break;
		int j = i + 1;
		while (j < ssize && cov_stats[j] >= cov) ++j;
		HBN_LOG("%d --- %d", i, j);
		i = j;
	}
	free(cov_stats);
}

int main(int argc, char* argv[])
{
	hbn_assert(argc == 3);
	CSeqDB* sdb = seqdb_load(argv[1], INIT_SUBJECT_DB_TITLE, 0);
	for (int i = 0; i < sdb->dbinfo.num_seqs; ++i) {
		cov_stat_for_one_contig(argv[2], sdb, i);
	}
	CSeqDBFree(sdb);
}
#endif

#if 1

#define E 100

#define NOL 0
#define LOL 1
#define ROL 2

static int
is_correct_ovlp(const int qb, const int qe, const int qs, 
    const int sb, const int se, const int ss)
{
	if (qb <= E && qs - qe <= E) return NOL;
	if (sb <= E && ss - se <= E) return NOL;
    if (ss - se <= E && qb <= E) return LOL;
    if (qs - qe <= E && sb <= E) return ROL;
    return NOL;
}

int main(int argc, char* argv[])
{
	ifstream in(argv[1]);
	string line;
	vector<ncbi::CTempString> components;
	const ncbi::CTempString kDelim("\t");
	vector<string> hlist;

	string last_name;
	int lol = 0;
	int rol = 0;

	while(getline(in, line)) {
		components.clear();
		NStr::Split(line, kDelim, components);
		hbn_assert(components.size() == 12, "%zu", components.size());

		int qdir = NStr::StringToInt(components[4]);
		int qoff = NStr::StringToInt(components[5]);
		int qend = NStr::StringToInt(components[6]);
		int qsize = NStr::StringToInt(components[7]);
		int soff = NStr::StringToInt(components[9]);
		int send = NStr::StringToInt(components[10]);
		int ssize = NStr::StringToInt(components[11]);
		double perc_identity = NStr::StringToDouble(components[2]);

		int qb, qe, sb, se;
		if (qdir == FWD) {
			qb = qoff;
			qe = qend;
			sb = soff;
			se = send;
		} else {
			hbn_assert(qdir == REV);
			qb = qsize - qend;
			qe = qsize - qoff;
			sb = ssize - send;
			se = ssize - soff;
		}
		string name = components[0];
		if (name != last_name) {
			if (lol && rol) {
				for (auto& e : hlist) cout << e << endl;
			}
			hlist.clear();
			last_name = name;
			lol = 0;
			rol = 0;
		}

		int r = is_correct_ovlp(qb, qe, qsize, sb, se, ssize);
		if (r) {
			hlist.push_back(line);
			if (r == LOL) lol = 1;
			if (r == ROL) rol = 1;
		}
	}
	in.close();
}
#endif