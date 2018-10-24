//
// Created by Amber on 2018/4/3.
//

#include "Layout.h"

void Layout::initLayout(fstream& inFile, InputMode& input_mode, int& helio_type) {
	stringstream line_stream;
	string line, word;
	while (getline(inFile, line)) {
		line_stream.clear();
		line_stream.str(line);
		line_stream >> word;
		if (word == "end") {
			input_mode = Initial;
			break;
		}
		else if (word == "pos")
			line_stream >> layout_bound_pos.x() >> layout_bound_pos.y() >> layout_bound_pos.z();
		else if (word == "size")
			line_stream >> layout_size.x() >> layout_size.y() >> layout_size.z();
		else if (word == "inter")
			line_stream >> helio_interval.x() >> helio_interval.y() >> helio_interval.z();
		else if (word == "n")
			line_stream >> helio_num;
		else
			line_stream >> helio_type;
	}
	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;
}

void Layout::setHelioLayout(vector<Heliostat*> helios)
{
	layout_first_helio_center.y() = helios[0]->helio_pos.y();
	layout_row_col.x() = layout_size.z() / helio_interval.z();		//row
	layout_row_col.y() = layout_size.x() / helio_interval.x();		//col
	helio_layout.resize(layout_row_col.x(), vector<vector<Heliostat*>>(layout_row_col.y()));

	vector<int> row_col(2, 0);
	for (auto&helio : helios) {
		set<vector<int>> pos;
		for (auto&v : helio->vertex) {
			row_col[0] = (v.z() - layout_bound_pos.z()) / helio_interval.z();	// smaller z is, smaller row is
			row_col[1] = (v.x() - layout_bound_pos.x()) / helio_interval.x();	// smaller x is, smaller col is
			// int i = pos.count(row_col);
			//pos.find(row_col);
			if (pos.count(row_col) == 0) {
				pos.insert(row_col);
				helio_layout[row_col[0]][row_col[1]].push_back(helio);
			}
		}
		
	}
}
