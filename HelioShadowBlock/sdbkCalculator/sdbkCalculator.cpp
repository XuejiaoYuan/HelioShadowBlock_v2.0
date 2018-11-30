#include"sdbkCalculator.h"


void SdBkCalc::sample_calc_preprocess(const int sample_row_num, const int sample_col_num, bool calc_s, bool calc_f)
{
	if (calc_f || field_index == nullptr)
		field_index = field_data_pre();
	if (calc_s || sample_field_index == nullptr)
		sample_field_index = sample_field_data_pre(sample_row_num, sample_col_num);
}

MatrixXf* SdBkCalc::calcSampleShadowBlock()
{
	vector<Receiver*> recvs = solar_scene->recvs;
	int row = sample_field_index->rows();
	int col = sample_field_index->cols();
	// Vector3f focus_center = recvs[0]->recv_pos + Vector3f(recvs[0]->recv_normal.array() * recvs[0]->recv_size.array());
	
	Vector3f reverse_sunray_dir = -solar_scene->sunray_dir;
	sample_sd_bk_res = new MatrixXf(row, col);

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			int index = (*sample_field_index)(i, j);
			Heliostat* helio = solar_scene->helios[index];
			set<vector<int>> shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda;
			// Vector3f reflect_dir = focus_center - helio->helio_pos;
			int fc_index = solar_scene->helios[index]->focus_center_index;
			Vector3f reflect_dir = recvs[0]->focus_center[fc_index] - helio->helio_pos;
			calcIntersection3DDDA(helio, reverse_sunray_dir, shadow_relative_grid_label_3ddda);

			calcIntersection3DDDA(helio, reflect_dir, block_relative_grid_label_3ddda);
			vector<Vector3f> dir = { reverse_sunray_dir, reflect_dir };
			vector<set<vector<int>>> estimate_grids = { shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda };
			(*sample_sd_bk_res)(i, j) = helioClipper(helio, dir, estimate_grids);
		}
	}
	return sample_sd_bk_res;
}


float SdBkCalc::helioClipper(Heliostat*helio, const Vector3f&dir, const set<vector<int>>& estimate_grid)
{
	vector<Vector3f> helio_v, local_v, tmp_v(4);
	vector<Vector2f> project_v(4);
	float t;
	Paths subj(1), clips;
	helio_v = helio->vertex;
	Vector3f reverse_dir = Vector3f(-dir.x(), -dir.y(), -dir.z());

	for (int i = 0; i < helio_v.size(); i++) {
		local_v.push_back(GeometryFunc::mulMatrix(helio_v[i], helio->world2localM));
		subj[0] << IntPoint(VERTEXSCALE*local_v[i].x(), VERTEXSCALE*local_v[i].z());
	}

	float total_area = 0;
	int subj_v_num = subj[0].size();
	for (int j = 0; j < subj_v_num; j++)
		total_area += (subj[0][j].X*subj[0][(j + 1) % subj_v_num].Y)
		- (subj[0][(j + 1) % subj_v_num].X*subj[0][j].Y);
	total_area = fabs(total_area*0.5);

	if (total_area == 0) {
		cout << "Project surface is 0!" << endl;
		return 0;
	}

	vector<Vector3f> pro(4);
	for (auto iter = estimate_grid.begin(); iter != estimate_grid.end(); iter++) {
		for (auto&relative_helio : solar_scene->layouts[0]->helio_layout[(*iter)[0]][(*iter)[1]]) {
			if (relative_helio == helio)
				continue;
			helio_v = relative_helio->vertex;
			int cnt = 0;
			for (int i = 0; i < helio_v.size(); i++) {
				t = calcIntersectionPoint(helio_v[i], reverse_dir, helio->vertex[0], helio->vertex[1], helio->vertex[2]);
				if (t < Epsilon)
					break;
				pro[i].x() = helio_v[i].x() + t*reverse_dir.x();
				pro[i].y() = helio_v[i].y() + t*reverse_dir.y();
				pro[i].z() = helio_v[i].z() + t*reverse_dir.z();
			}
			if (t >= Epsilon) {
				Path clip;
				for (auto v : pro) {
					v = GeometryFunc::mulMatrix(v, helio->world2localM);
					clip << IntPoint(VERTEXSCALE *v.x(), VERTEXSCALE * v.z());
				}
				clips.push_back(clip);
			}
		}
	}

	Clipper c;
	Paths solution;													// solution represents the shadowing / blocking area
	c.AddPaths(subj, ptSubject, true);
	c.AddPaths(clips, ptClip, true);
	c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);

	float sum = 0;
	for (int i = 0; i < solution.size(); i++) {
		int n = solution[i].size();
		for (int j = 0; j < n; j++) {
			sum += (solution[i][j].X*solution[i][(j + 1) % n].Y)
				- (solution[i][(j + 1) % n].X*solution[i][j].Y);
		}
	}
	sum = fabs(sum*0.5);

	float res = sum / total_area;

#ifdef OUTPUTRES
	fstream outFile;
	outFile.open("python_clipper.txt", ios_base::out | ios_base::app);
	outFile << "subj:" << endl;
	for (auto&tmp : subj[0])
		outFile << tmp.X << ' ' << tmp.Y << endl;
	outFile << "clips:" << endl;
	for (int i = 0; i < clips.size(); i++)
		for (int j = 0; j < clips[i].size(); j++)
			outFile << clips[i][j].X << ' ' << clips[i][j].Y << endl;
	outFile << "solution:" << endl;
	for (int i = 0; i < solution.size(); i++)
		for (int j = 0; j < solution[i].size(); j++)
			outFile << solution[i][j].X << ' ' << solution[i][j].Y << endl;
	outFile.close();

#endif // OUTPUTRES

	return res;
}

float SdBkCalc::helioClipper(Heliostat * helio, const vector<Vector3f>& dir, const vector<set<vector<int>>>& estimate_grids)
{
	vector<Vector3f> helio_v, local_v, tmp_v(4);
	vector<Vector2f> project_v(4);
	float t;
	Paths subj(1), clips;
	helio_v = helio->vertex;

	for (int i = 0; i < helio_v.size(); i++) {
		local_v.push_back(GeometryFunc::mulMatrix(helio_v[i], helio->world2localM));
		subj[0] << IntPoint(VERTEXSCALE*local_v[i].x(), VERTEXSCALE*local_v[i].z());
	}

	long long total_area = 0;
	int subj_v_num = subj[0].size();
	for (int j = 0; j < subj_v_num; j++)
		total_area += (subj[0][j].X*subj[0][(j + 1) % subj_v_num].Y)
		- (subj[0][(j + 1) % subj_v_num].X*subj[0][j].Y);
	total_area = fabs(total_area*0.5);

	if (total_area == 0) {
		cout << "Project surface is 0!" << endl;
		return 0;
	}

	for (int index = 0; index < 2; index++) {
		Vector3f reverse_dir = Vector3f(-dir[index].x(), -dir[index].y(), -dir[index].z());
		vector<Vector3f> pro(4);
		set<Heliostat*> relative_helio_set;
		for (auto iter = estimate_grids[index].begin(); iter != estimate_grids[index].end(); iter++) {
			for (auto&relative_helio : solar_scene->layouts[0]->helio_layout[(*iter)[0]][(*iter)[1]]) {
				if (relative_helio_set.count(relative_helio) || relative_helio == helio)
					continue;
				else
					relative_helio_set.insert(relative_helio);
				helio_v = relative_helio->vertex;
				int cnt = 0;
				for (int i = 0; i < helio_v.size(); i++) {
					t = calcIntersectionPoint(helio_v[i], reverse_dir, helio->vertex[0], helio->vertex[1], helio->vertex[2]);
					if (t < Epsilon)
						break;
					pro[i].x() = helio_v[i].x() + t*reverse_dir.x();
					pro[i].y() = helio_v[i].y() + t*reverse_dir.y();
					pro[i].z() = helio_v[i].z() + t*reverse_dir.z();
				}
				if (t >= Epsilon) {
					Path clip;
					for (auto v : pro) {
						v = GeometryFunc::mulMatrix(v, helio->world2localM);
						clip << IntPoint(VERTEXSCALE *v.x(), VERTEXSCALE * v.z());
					}
					clips.push_back(clip);
				}
			}
		}
	}

	Clipper c;
	Paths solution;													// solution represents the shadowing / blocking area
	c.AddPaths(subj, ptSubject, true);
	c.AddPaths(clips, ptClip, true);
	c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);

	float sum = 0;
	for (int i = 0; i < solution.size(); i++) {
		int n = solution[i].size();
		for (int j = 0; j < n; j++) {
			sum += (solution[i][j].X*solution[i][(j + 1) % n].Y)
				- (solution[i][(j + 1) % n].X*solution[i][j].Y);
		}
	}
	sum = fabs(sum*0.5);

	float res = sum / total_area;

#ifdef OUTPUTRES
	fstream outFile;
	outFile.open("python_clipper.txt", ios_base::out | ios_base::app);
	outFile << "subj:" << endl;
	for (auto&tmp : subj[0])
		outFile << tmp.X << ' ' << tmp.Y << endl;
	outFile << "clips:" << endl;
	for (int i = 0; i < clips.size(); i++)
		for (int j = 0; j < clips[i].size(); j++)
			outFile << clips[i][j].X << ' ' << clips[i][j].Y << endl;
	outFile << "solution:" << endl;
	for (int i = 0; i < solution.size(); i++)
		for (int j = 0; j < solution[i].size(); j++)
			outFile << solution[i][j].X << ' ' << solution[i][j].Y << endl;
	outFile.close();

#endif // OUTPUTRES

	return res;
}


//	使用ray tracing方式计算blocking与shadowing
//	将定日镜均匀剖分，每个小区域中心点 发出一根光线，与周围定日镜求交
//	记录产生blocking和shadowing的定日镜的编号
//	version: CPU
float SdBkCalc::calcAccurateIntersection(Heliostat * helio, const Vector3f & dir, set<vector<int>>& relative_helio_label)
{
	int helio_inner_rows = 20;
	int helio_inner_cols = 20;

	float total_sum = 0;
	int cnt = 0;

	vector<Vector3f> subhelio_v;
	helio->getSubHelioVertex(subhelio_v);
	int subHelio_num = subhelio_v.size() / 4;
	for (int i = 0; i < subHelio_num; i++) {
		vector<Vector3f> helio_v;
		for (int k = 0; k < 4; k++)
			helio_v.push_back(subhelio_v[4 * i + k]);
		Vector3f row_dir = (helio_v[1] - helio_v[0]) / helio_inner_rows;
		Vector3f col_dir = (helio_v[3] - helio_v[0]) / helio_inner_cols;

		for (int i = 0; i <= helio_inner_rows; i++) {
			for (int j = 0; j <= helio_inner_cols; j++) {
				total_sum++;
				Vector3f ori_v = helio_v[0] + i*row_dir + j*col_dir;
				float tMin = INT_MAX;
				Heliostat* hNear = nullptr;
				for (auto&h : solar_scene->helios) {
					if (h != helio) {
						vector<Vector3f> neigh_v = h->vertex;

						Vector3f E1 = neigh_v[1] - neigh_v[0];
						Vector3f E2 = neigh_v[2] - neigh_v[0];
						Vector3f pvec = dir.cross(E2);
						float det = E1.dot(pvec);

						// ray and triangle are parallel if det is close to 0
						if (fabsf(det) < Epsilon) continue;

						float invDet = 1 / det;

						Vector3f T = ori_v - neigh_v[0];
						Vector3f qvec = T.cross(E1);

						float t = E2.dot(qvec)*invDet;
						if (t < Epsilon) continue;

						Vector3f intersect_v = ori_v + dir*t;
						int l;
						for (l = 0; l < 4; l++) {
							Vector3f edg = neigh_v[(l + 1) % 4] - neigh_v[l];
							Vector3f line = intersect_v - neigh_v[l];
							Vector3f tmp_n = edg.cross(line);
							if (tmp_n.dot(Vector3f(0, 1, 0)) < Epsilon)
								break;
						}
						if (l != 4)
							continue;

						if (t < tMin) {
							tMin = t;
							hNear = h;
						}
					}
				}
				if (hNear != nullptr) {
					//vector<int> res = { (int)hNear->helio_label.x(), (int)hNear->helio_label.y(), (int)hNear->helio_label.z() };
					vector<int> res(3);
					res[0] = (hNear->helio_pos.z() - solar_scene->layouts[0]->layout_first_helio_center.z()) / solar_scene->layouts[0]->helio_interval.z();		// smaller x is, smaller col is
					res[1] = (hNear->helio_pos.x() - solar_scene->layouts[0]->layout_first_helio_center.x()) / solar_scene->layouts[0]->helio_interval.x();		// bigger z is, smaller row is
					res[2] = 0;
					relative_helio_label.insert(res);
					cnt++;
				}
			}
		}
	}
	float res = cnt / total_sum;

	return res;
}


float SdBkCalc::calcAccurateIntersection(Heliostat* helio, const vector<Vector3f>& dir, vector<set<vector<int>>>& relative_helio_label)
{
	int helio_inner_rows = 20;
	int helio_inner_cols = 20;

	float total_sum = 0;
	int cnt = 0;

	//vector<Vector3f> subhelio_v;
	//helio->getSubHelioVertex(subhelio_v);
	//int subHelio_num = subhelio_v.size() / 4;
	for (int i = 0; i < helio->vertex.size(); i++) {
		vector<Vector3f> helio_v = helio->vertex;
		//vector<Vector3f> helio_v;
		//for (int k = 0; k < 4; k++)
		//	helio_v.push_back(subhelio_v[4 * i + k]);
		Vector3f row_dir = (helio_v[1] - helio_v[0]) / helio_inner_rows;
		Vector3f col_dir = (helio_v[3] - helio_v[0]) / helio_inner_cols;

		for (int i = 0; i <= helio_inner_rows; i++) {
			for (int j = 0; j <= helio_inner_cols; j++) {
				total_sum++;
				Vector3f ori_v = helio_v[0] + i*row_dir + j*col_dir;
				float tMin = INT_MAX;
				Heliostat* hNear = nullptr;
				for (auto&h : solar_scene->helios) {
					if (h != helio) {
						vector<Vector3f> neigh_v = h->vertex;

						Vector3f E1 = neigh_v[1] - neigh_v[0];
						Vector3f E2 = neigh_v[2] - neigh_v[0];
						Vector3f pvec = dir[0].cross(E2);
						float det = E1.dot(pvec);

						// ray and triangle are parallel if det is close to 0
						if (fabsf(det) < Epsilon) continue;

						float invDet = 1 / det;

						Vector3f T = ori_v - neigh_v[0];
						Vector3f qvec = T.cross(E1);

						float t = E2.dot(qvec)*invDet;
						if (t < Epsilon) continue;

						Vector3f intersect_v = ori_v + dir[0]*t;
						int l;
						for (l = 0; l < 4; l++) {
							Vector3f edg = neigh_v[(l + 1) % 4] - neigh_v[l];
							Vector3f line = intersect_v - neigh_v[l];
							Vector3f tmp_n = edg.cross(line);
							if (tmp_n.dot(Vector3f(0, 1, 0)) < Epsilon)
								break;
						}
						if (l != 4)
							continue;

						if (t < tMin) {
							tMin = t;
							hNear = h;
						}
					}
				}
				if (hNear != nullptr) {
					vector<int> res(3);
					res[0] = (hNear->helio_pos.z() - solar_scene->layouts[0]->layout_first_helio_center.z()) / solar_scene->layouts[0]->helio_interval.z();		// smaller x is, smaller col is
					res[1] = (hNear->helio_pos.x() - solar_scene->layouts[0]->layout_first_helio_center.x()) / solar_scene->layouts[0]->helio_interval.x();		// bigger z is, smaller row is
					res[2] = 0;
					relative_helio_label[0].insert(res);
					cnt++;
				}
				else {
					for (auto&h : solar_scene->helios) {
						if (h != helio) {
							vector<Vector3f> neigh_v = h->vertex;

							Vector3f E1 = neigh_v[1] - neigh_v[0];
							Vector3f E2 = neigh_v[2] - neigh_v[0];
							Vector3f pvec = dir[1].cross(E2);
							float det = E1.dot(pvec);

							// ray and triangle are parallel if det is close to 0
							if (fabsf(det) < Epsilon) continue;

							float invDet = 1 / det;

							Vector3f T = ori_v - neigh_v[0];
							Vector3f qvec = T.cross(E1);

							float t = E2.dot(qvec)*invDet;
							if (t < Epsilon) continue;

							Vector3f intersect_v = ori_v + dir[1] * t;
							int l;
							for (l = 0; l < 4; l++) {
								Vector3f edg = neigh_v[(l + 1) % 4] - neigh_v[l];
								Vector3f line = intersect_v - neigh_v[l];
								Vector3f tmp_n = edg.cross(line);
								if (tmp_n.dot(Vector3f(0, 1, 0)) < Epsilon)
									break;
							}
							if (l != 4)
								continue;

							if (t < tMin) {
								tMin = t;
								hNear = h;
							}
						}
					}
					if (hNear != nullptr) {
						vector<int> res(3);
						res[0] = (hNear->helio_pos.z() - solar_scene->layouts[0]->layout_first_helio_center.z()) / solar_scene->layouts[0]->helio_interval.z();		// smaller x is, smaller col is
						res[1] = (hNear->helio_pos.x() - solar_scene->layouts[0]->layout_first_helio_center.x()) / solar_scene->layouts[0]->helio_interval.x();		// bigger z is, smaller row is
						res[2] = 0;
						relative_helio_label[1].insert(res);

						cnt++;
					}
				}
			}
		}
	}
	float res = cnt / total_sum;

	return res;
}


// 使用3DDDA计算产生阴影和遮挡的相关定日镜
// version: CPU
void SdBkCalc::calcIntersection3DDDA(Heliostat * helio, const Vector3f & dir, set<vector<int>>& relative_helio_label)
{
	vector<Layout*>& layouts = solar_scene->layouts;
	vector<Vector3f> helio_v;
	vector<Vector2f> project_helio_v;

	helio_v = helio->vertex;
	for (int i = 0; i < 4; i++) {
		project_helio_v.push_back(Vector2f(helio_v[i].x(), helio_v[i].z()));
	}

	// 1. 確定光線在場地中y方向移動距離及最終離開網格的位置
	float upper_y = layouts[0]->layout_size.y() + layouts[0]->layout_bound_pos.y();
	Vector2f upper_v[4];
	for (int i = 0; i < 4; i++) {
		float dis = (upper_y - helio_v[i].y()) / dir.y();
		upper_v[i] = Vector2f(dis*dir.x(), dis*dir.z()) + project_helio_v[i];
	}

	// 2. 確定定日鏡反射光柱在鏡場中的矩形包圍盒範圍
	vector<Vector2f> boundBox(2);
	boundBox[0] = Vector2f(INT_MAX, INT_MAX);		// min boundary
	boundBox[1] = Vector2f(INT_MIN, INT_MIN);		// max boundary
	for (int i = 0; i < 4; i++) {
		boundBox[0].x() = fmin(boundBox[0].x(), fmin(project_helio_v[i].x(), upper_v[i].x()));
		boundBox[0].y() = fmin(boundBox[0].y(), fmin(project_helio_v[i].y(), upper_v[i].y()));
		boundBox[1].x() = fmax(boundBox[1].x(), fmax(project_helio_v[i].x(), upper_v[i].x()));
		boundBox[1].y() = fmax(boundBox[1].y(), fmax(project_helio_v[i].y(), upper_v[i].y()));
	}
	Vector2f ray_dir = Vector2f(dir.x(), dir.z()).normalized();
	Vector2f deltaT;
	Vector2f t;
	Vector2f origs[2] = {
		Vector2f(INT_MAX,INT_MAX),		// min vertex
		Vector2f(INT_MIN,INT_MIN)		// max vertex
	};

	
	// 3. 確定layout下各網格的矩陣間隔
	Vector2f cellDimension = Vector2f(layouts[0]->helio_interval.x(), layouts[0]->helio_interval.z());

	float colLength = boundBox[1].x() - boundBox[0].x();
	float rowLength = boundBox[1].y() - boundBox[0].y();

	int helio_col = (helio->helio_pos.x() - layouts[0]->layout_first_helio_center.x()) / cellDimension.x();			// smaller x is, smaller col is
	int helio_row = (helio->helio_pos.z() - layouts[0]->layout_first_helio_center.z()) / cellDimension.y();			// smaller z is, smaller row is

	int minCol = (boundBox[0].x() - layouts[0]->layout_first_helio_center.x()) / cellDimension.x();
	int minRow = (boundBox[0].y() - layouts[0]->layout_first_helio_center.z()) / cellDimension.y();
	int maxCol = (boundBox[1].x() - layouts[0]->layout_first_helio_center.x()) / cellDimension.x() + 0.5;
	int maxRow = (boundBox[1].y() - layouts[0]->layout_first_helio_center.z()) / cellDimension.y() + 0.5;

	minCol = max(0, minCol);
	minRow = max(0, minRow);
	maxCol = min(layouts[0]->layout_row_col.y() - 1, maxCol);
	maxRow = min(layouts[0]->layout_row_col.x() - 1, maxRow);

	for (auto&orig : project_helio_v) {
		Vector2f o_grid(
			orig.x() - layouts[0]->layout_bound_pos.x(),
			orig.y() - layouts[0]->layout_bound_pos.z()
		);

		if (ray_dir.x() < 0) {
			deltaT.x() = -cellDimension.x() / ray_dir.x();
			t.x() = (floor(o_grid.x() / cellDimension.x())*cellDimension.x() - o_grid.x()) / ray_dir.x();
		}
		else {
			deltaT.x() = cellDimension.x() / ray_dir.x();
			t.x() = ((floor(o_grid.x() / cellDimension.x()) + 1)*cellDimension.x() - o_grid.x()) / ray_dir.x();
		}

		if (ray_dir.y() < 0) {
			deltaT.y() = -cellDimension.y() / ray_dir.y();
			t.y() = (floor(o_grid.y() / cellDimension.y())*cellDimension.y() - o_grid.y()) / ray_dir.y();
		}
		else {
			deltaT.y() = cellDimension.y() / ray_dir.y();
			t.y() = ((floor(o_grid.y() / cellDimension.y()) + 1)*cellDimension.y() - o_grid.y()) / ray_dir.y();
		}

		float tmp = 0;
		Vector2i grid_label(
			(orig.y() - layouts[0]->layout_bound_pos.z()) / cellDimension.y(),	// smaller z is, smaller row is
			(orig.x() - layouts[0]->layout_bound_pos.x()) / cellDimension.x()	// smaller x is, smaller col is
		);

		while (1) {
			if (grid_label.x() < minRow || grid_label.x() > maxRow ||
				grid_label.y() < minCol || grid_label.y() > maxCol)
				break;
			else if (grid_label.x() != helio_row || grid_label.y() != helio_col) {
				vector<int> res = { (int)grid_label.x(), (int)grid_label.y(), 0 };
				relative_helio_label.insert(res);
			}

			if (t.x() < t.y()) {
				tmp = t.x();
				t.x() += deltaT.x();
				if (ray_dir.x() < 0)
					grid_label.y()--;
				else
					grid_label.y()++;
			}
			else {
				tmp = t.y();
				t.y() += deltaT.y();
				if (ray_dir.y() < 0)		// smaller z is, smaller row is
					grid_label.x()--;
				else
					grid_label.x()++;
			}
		}
	}
}


float SdBkCalc::calcIntersectionPoint(const Vector3f & orig, const Vector3f & dir, const Vector3f & A, const Vector3f & B, const Vector3f & C)
{
	Vector3f E1 = B - A;
	Vector3f E2 = C - A;
	Vector3f pvec = dir.cross(E2);
	float det = E1.dot(pvec);

	// ray and triangle are parallel if det is close to 0
	if (fabsf(det) < Epsilon) return 0;

	float invDet = 1 / det;

	Vector3f T = orig - A;

	Vector3f qvec = T.cross(E1);

	float t = E2.dot(qvec)*invDet;
	return t;
}


 
// helio: 待计算定日镜
// dir: 定日镜到接收器的反射光线方向
float SdBkCalc::calcFluxMap(Heliostat * helio, const float DNI)
{
	//auto start = std::chrono::high_resolution_clock::now();


	int fc_index = helio->focus_center_index;
	vector<Receiver*> recvs = solar_scene->recvs;
	Vector3f focus_center = recvs[0]->focus_center[fc_index];
	Vector3f reverse_dir = (helio->helio_pos - focus_center).normalized();		// The normal of image plane
	float flux_sum = 0;

	Matrix4f word2localM, local2worldM;
	GeometryFunc::getMatrixs(reverse_dir, focus_center, local2worldM, word2localM);

	for (int i = 0; i < helio->cos_phi.size(); i++) {
		if (helio->cos_phi[i] > 0) {
			vector<Vector2d> proj_v;
			for (auto& v : recvs[0]->recv_vertex[i]) {
				Vector3f inter_v = GeometryFunc::calcIntersection(reverse_dir, focus_center, v, reverse_dir);
				//cout << "pos: " << v.x() << ' ' << v.y()<< ' ' << v.z() << endl;
				//cout << "inter: " << inter_v.x() << ' ' <<inter_v.y() << ' '  << inter_v.z() << endl;
				
				inter_v = GeometryFunc::mulMatrix(inter_v, word2localM);
				proj_v.push_back(Vector2d(inter_v.x(), inter_v.z()));
				//cout << "proj: " << inter_v.x() << ' ' << inter_v.z() << endl;

			}
			// flux_sum += _calc_flux_sum(proj_v, recvs[0]->mask_rows, recvs[0]->mask_cols, helio, helio->cos_phi[i], DNI);
			float _flux_sum_grid1 = _calc_flux_sum(proj_v, recvs[0]->mask_rows, recvs[0]->mask_cols, helio, helio->cos_phi[i], DNI);
			float _flux_sum_grid2 = _calc_flux_sum(proj_v, helio, helio->cos_phi[i], DNI);
			// cout << _flux_sum_grid1 << ' ' << _flux_sum_grid2 << endl;
		}
	}
	//auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	//auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	//std::cout << "Calculate flux sum time: " << time << "s." << endl;
	
	return flux_sum;
}

///
// 对接收器投影面进行均匀分割，统计每个分割面的flux分布结果，并求和
///
inline float SdBkCalc::_calc_flux_sum(vector<Vector2d>& proj_v,const int rows, const int cols, Heliostat * helio, const float cos_phi, const float DNI)
{

	auto start = std::chrono::high_resolution_clock::now();

	MatrixXf mask_x(rows, cols), mask_y(rows, cols);
	float row_gap = (proj_v[2] - proj_v[1]).y() / rows;			// 矩阵列方向变化率
	float col_gap = (proj_v[0] - proj_v[1]).x() / cols;			// 矩阵行方向变化率
	for (int i = 0; i < rows; i++) {
		float start_v = proj_v[1].x() + i* col_gap;
		float end_v = proj_v[2].x() + i*col_gap;
		mask_x.row(i).setLinSpaced(cols, start_v, end_v);
	}

	for (int i = 0; i < cols; i++) {
		float start_v = proj_v[1].y() + i * row_gap;
		float end_v = proj_v[0].y() + i* row_gap;
		mask_y.col(i).setLinSpaced(rows, 1/ helio->l_w_ratio * start_v, 1/ helio->l_w_ratio * end_v);
	}

	float flux_sum = (-0.5 / pow(helio->sigma, 2) * (mask_x.array().pow(2) + mask_y.array().pow(2))).array().exp().sum();
	cout << "grid: " << flux_sum << endl;
	flux_sum = flux_sum * helio->flux_param * (1 - helio->sd_bk) * DNI * cos_phi;

	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	// std::cout << "Calculate flux sum time: " << time << "s." << endl;

	return flux_sum;
}

///
// 对接收器投影面进行Gauss-Legendre数值积分
///
inline float SdBkCalc::_calc_flux_sum(vector<Vector2d>& proj_v, Heliostat * helio, const float cos_phi, const float DNI)
{
	auto start = std::chrono::high_resolution_clock::now();

	vector<VectorXd> weight = gl->getW();
	vector<VectorXd> node = gl->getX();
	double k1 = INFINITY;
	double k2 = INFINITY;
	double u1, u2, v1, v2;
	if (abs(proj_v[1].x() - proj_v[0].x()) < Epsilon) {
		u1 = proj_v[0].x();
		u2 = proj_v[3].x();
	}
	else {
		k1 = (proj_v[1].y() - proj_v[0].y()) / (proj_v[1].x() - proj_v[0].x());
		u1 = k1*proj_v[0].x() + proj_v[0].y();
		u2 = k1*proj_v[3].x() + proj_v[3].y();
	}
	if (abs(proj_v[2].x() - proj_v[1].x()) < Epsilon) {
		v1 = proj_v[0].x();
		v2 = proj_v[1].x();
	}
	else {
		k2 = (proj_v[2].y() - proj_v[1].y()) / (proj_v[2].x() - proj_v[1].x());
		v1 = k2*proj_v[0].x() + proj_v[0].y();
		v2 = k2*proj_v[1].x() + proj_v[1].y();
	}

	double a = min(u1, u2);
	double b = max(u1, u2);
	double c = min(v1, v2);
	double d = max(v1, v2);
	double alpha = (d - c) / 2;
	double beta = (b - a) / 2;
	double item1 = (b + a) / 2;
	double item2 = (d + c) / 2;

	double sum = 0;
	double u, v;
	double mute = 1.0;
	for (int i = 0; i < weight[0].size(); i++)
		cout << weight[0](i) << endl;
	cout << endl;
	for (int i = 0; i < weight[1].size(); i++)
		cout << weight[1](i) << endl;
	for (int i = 0; i < node[0].size(); i++) {
		for (int j = 0; j < node[1].size(); j++) {
			u = node[0](i)*beta + item1;
			v = node[1](j)*alpha + item2;
			if (k1 == INFINITY && abs(k2) > Epsilon) {
				v = v - k2*u;
			}
			else {
				mute = k1 - k2;
				v = (u - v) / (k1 - k2);
				u = (k1*v - k2*u) / (k1 - k2);
			}
			// cout << v << ' ' << u << endl;
			sum+= mute * alpha* beta* weight[0](i)*weight[1](j)
				*exp(-0.5 / pow(helio->sigma, 2)*(pow(u, 2) + pow(v/helio->l_w_ratio, 2)));
		}
	}
	cout << "integer: " << sum << endl;
	sum = sum * helio->flux_param * (1 - helio->sd_bk) * DNI * cos_phi;
	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	// std::cout << "Calculate flux sum time: " << time << "s." << endl;


	return sum;
}


// 以光线跟踪计算结果为groundtruth，检测预测结果是否包含准确结果
// version: CPU
float SdBkCalc::checkForRelativeHelio(const set<vector<int>>& accurate_helio, const set<vector<int>>& estimate_helio)
{
	// check for estimation of relative heliostat
	int sum = estimate_helio.size();
	int cnt = accurate_helio.size();
	bool flag = true;
	for (int helio_row = 0; helio_row < solar_scene->layouts[0]->layout_row_col.x(); helio_row++)
		for (int helio_col = 0; helio_col < solar_scene->layouts[0]->layout_row_col.y(); helio_col++) {
			for (auto it = accurate_helio.begin();
				it != accurate_helio.end(); it++) {
				if (!estimate_helio.count(*it))
				{
					cnt--;
					cout << "row: " << helio_row << " col: " << helio_col << endl;
					cout << "Didn't contain row: " << (*it)[0] << " col: " << (*it)[1] << endl;
				}
			}
		}
	//return sum > 0 ? (float)cnt / (float)sum : 0;
	return (float)cnt / (float)solar_scene->helios.size();
}


float SdBkCalc::calcSingleShadowBlock(int helio_index)
{
	vector<Receiver*> recvs = solar_scene->recvs;

	auto helio = this->solar_scene->helios[helio_index];
	set<vector<int>> shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda;
	Vector3f reverse_sunray_dir = -solar_scene->sunray_dir;
	// Vector3f focus_center = recvs[0]->recv_pos + Vector3f(recvs[0]->recv_normal.array() * recvs[0]->recv_size.array());
	int fc_index = solar_scene->helios[helio_index]->focus_center_index;
	Vector3f reflect_dir = recvs[0]->focus_center[fc_index] - helio->helio_pos;

	// 3DDDA + clipper 
	//Calc the relative heliostats which cause shadowing
	calcIntersection3DDDA(helio, reverse_sunray_dir, shadow_relative_grid_label_3ddda);

	//Calc the relactive heliostats which cause blocking
	// Vector3f reflect_dir = (focus_center - helio->helio_pos).normalized();
	calcIntersection3DDDA(helio, reflect_dir, block_relative_grid_label_3ddda);

	vector<Vector3f> dir = { reverse_sunray_dir, reflect_dir };
	vector<set<vector<int>>> estimate_grids = { shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda };
	return helioClipper(helio, dir, estimate_grids);
}


void SdBkCalc::calcShadowBlock(const float DNI)
{
	vector<Heliostat*> helios = solar_scene->helios;
	vector<Layout*> layouts = solar_scene->layouts;
	vector<Receiver*> recvs = solar_scene->recvs;

	//int row = layouts[0]->layout_row_col.x();
	//int col = layouts[0]->layout_row_col.y();
	//sd_bk_res = new MatrixXf(row, col);

//#ifdef DEBUG
//	MatrixXf* raytracing_store = new MatrixXf(row, col);
//#endif // DEBUG


#ifdef READFILE
	fstream inFile("shadowblock_gt_save.txt", ios_base::in);
	if (inFile.fail())
		cout << "Can't open the file!" << endl;
#endif // READFILE

	auto start = std::chrono::high_resolution_clock::now();

	Vector3f reverse_sunray_dir = -solar_scene->sunray_dir;
#pragma omp parallel for
	for (int i = 0; i < helios.size(); i++) {

		auto helio = helios[i];
		set<vector<int>> shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda;
		//int r, c;
		//get_row_col(helio->helio_index, r, c);

		// 3DDDA + clipper 
		//Calc the relative heliostats which cause shadowing
		calcIntersection3DDDA(helio, reverse_sunray_dir, shadow_relative_grid_label_3ddda);

		//Calc the relactive heliostats which cause blocking
		// Vector3f reflect_dir = (focus_center - helio->helio_pos).normalized();
		int fc_index = helio->focus_center_index;
		Vector3f reflect_dir = (recvs[0]->focus_center[fc_index] - helio->helio_pos).normalized();

		calcIntersection3DDDA(helio, reflect_dir, block_relative_grid_label_3ddda);

		vector<Vector3f> dir = { reverse_sunray_dir, reflect_dir };
		if (shadow_relative_grid_label_3ddda.size() != 0 || block_relative_grid_label_3ddda.size() != 0) {
			vector<set<vector<int>>> estimate_grids = { shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda };
			helio->sd_bk = helioClipper(helio, dir, estimate_grids);
		}
		else
			helio->sd_bk = 0;
		
		if(gl!=NULL)
			helio->flux_sum = calcFluxMap(helio, DNI);

#ifdef DEBUG
		// Ray tracing
		vector<set<vector<int>>> ray_relative_grids(2);
		int helio_col = (helio->helio_pos.x() - layouts[0]->layout_first_helio_center.x()) / layouts[0]->helio_interval.x();			// smaller x is, smaller col is
		int helio_row = (helio->helio_pos.z() - layouts[0]->layout_first_helio_center.z()) / layouts[0]->helio_interval.z();

		std::cout << "Heliostat: " << helio->helio_index << endl;

		float ray_res = calcAccurateIntersection(helio, dir, ray_relative_grids);
		cout << ray_res - helio->sd_bk << endl;
		int r1 = (helio->helio_pos.z() - solar_scene->layouts[0]->layout_first_helio_center.z()) / solar_scene->layouts[0]->helio_interval.z();		// smaller x is, smaller col is
		int c1 = (helio->helio_pos.x() - solar_scene->layouts[0]->layout_first_helio_center.x()) / solar_scene->layouts[0]->helio_interval.x();		// bigger z is, smaller row is
		cout << r1 << ' ' << c1 << endl;
		//(*raytracing_store)(r,c) = calcAccurateIntersection(helio, dir);
		//cout << (*raytracing_store)(r, c) - (*sd_bk_res)(r, c) << endl;
		//for (auto&label : tmp_block_relative_grid_label) {
		//	if (block_relative_grid_label_3ddda.count(label) == 0) {
		//		cout << helio_col << ' ' << helio_row << ' ' << helio->helio_index << endl;
		//	}
		//}
#endif // DEBUG

#ifdef READFILE
		if (helio->helio_label.y() == 0)
			std::cout << "Heliostat: " << helio->helio_label.x() << endl;
		inFile >> raytracing_store[i][0] >> raytracing_store[i][1];
		cout << raytracing_store[i][0] - shadow_helio_ratio << ' ' << raytracing_store[i][1] - block_helio_ratio << endl;
		gt_field_shadow_ratio += raytracing_store[i][0];
		gt_field_block_ratio += raytracing_store[i][1];
#endif // READFILE
	}

	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	std::cout << "Heliostat clipper time: " << time << "s." << endl;


#ifdef OUTPUTRES
	fstream outFile("shadowblock_clipper_save.txt", ios_base::out);
	for (auto&res : clipper_res_store)
		outFile << res[0] << ' ' << res[1] << endl;
	outFile.close();

#endif // OUTPUTRES

	//return sd_bk_res;
}

void CrossRectSdBkCalc::save_clipper_res(const string save_path, int month, int day, int hour, int minute)
{
	//fstream outFile(save_path + "clipper_m" + to_string(month) + "_d" + to_string(day) + "_h" + to_string(hour) + "_min" + to_string(minute) + ".txt", ios_base::out);
	//int row = sd_bk_res->rows();
	//int col = sd_bk_res->cols();
	////int row = clipper_res_store[0]->rows();
	////int col = clipper_res_store[0]->cols();
	//int tmp_col;
	//outFile << row << ' ' << col << endl;
	//int cnt = 0;
	//for (int i = 0; i < row; i++) {
	//	if (i % 2) tmp_col = col - 1;
	//	else tmp_col = col;
	//	for (int j = 0; j < tmp_col; j++) {
	//		outFile << solar_scene->helios[cnt]->helio_pos.x() << ' '
	//			<< solar_scene->helios[cnt]->helio_pos.z() << ' '
	//			<< (*sd_bk_res)(i, j) << endl;
	//			//<< (*clipper_res_store[0])(i, j) << ' '
	//			//<< (*clipper_res_store[1])(i, j) << endl;
	//		cnt++;
	//	}
	//}
	//outFile.close();
}

void CrossRectSdBkCalc::get_row_col(const int index, int & r, int & c)
{
	int col = solar_scene->layouts[0]->layout_row_col.y();
	int g_cnts = 2 * col - 1;
	r = 2 * (index / g_cnts) + (index%g_cnts) / col;
	c = index%g_cnts%col;
}


MatrixXf* CrossRectSdBkCalc::field_data_pre() {
	vector<Heliostat*>& helios = solar_scene->helios;
	vector<Layout*>& layouts = solar_scene->layouts;
	int row = layouts[0]->layout_row_col.x();
	int col = layouts[0]->layout_row_col.y();
	if (!field_data.empty()) {
		delete field_data[0];
		delete field_data[1];
		field_data.clear();
	}

	MatrixXf* field_data_x = new MatrixXf(row, col);
	MatrixXf* field_data_y = new MatrixXf(row, col);
	MatrixXf* field_index = new MatrixXf(row, col);
	int h = 0;
	int tmp_col;
	for (int i = 0; i < row; i++) {
		if (i % 2) tmp_col = col - 1;
		else tmp_col = col;
		for (int j = 0; j < tmp_col; j++) {
			(*field_data_x)(i, j) = helios[h]->helio_pos.x();
			(*field_data_y)(i, j) = helios[h]->helio_pos.z();
			(*field_index)(i, j) = helios[h]->helio_index;
			h++;
		}
	}

	field_data.push_back(field_data_x);
	field_data.push_back(field_data_y);
	return field_index;
}


MatrixXf* CrossRectSdBkCalc::sample_field_data_pre(const int sample_row_num, const int sample_col_num)
{
	int row = solar_scene->layouts[0]->layout_row_col.x();
	int col = solar_scene->layouts[0]->layout_row_col.y();
	if (!sample_field_data.empty()) {
		delete sample_field_data[0];
		delete sample_field_data[1];
		sample_field_data.clear();
	}

	MatrixXf* sample_f_data_x = new MatrixXf(sample_row_num, sample_col_num);
	MatrixXf* sample_f_data_y = new MatrixXf(sample_row_num, sample_col_num);
	MatrixXf* sample_index = new MatrixXf(sample_row_num, sample_col_num);
	int tmp_col;
	int cur_row = 0;
	for (int i = 0; i < sample_row_num / 2 - 1; i++) {
		int f_i = int(row*i / sample_row_num);
		for (int k = 0; k < 2; k++) {
			for (int j = 0; j < sample_col_num - 1; j++) {
				int f_j = int(col*j / sample_col_num);
				(*sample_f_data_x)(cur_row, j) = (*field_data[0])(2 * f_i + k, f_j);
				(*sample_f_data_y)(cur_row, j) = (*field_data[1])(2 * f_i + k, f_j);
				(*sample_index)(cur_row, j) = (*field_index)(2 * f_i + k, f_j);
			}
			if (k % 2) tmp_col = col - 2;
			else tmp_col = col - 1;
			(*sample_f_data_x)(cur_row, sample_col_num - 1) = (*field_data[0])(2 * f_i + k, tmp_col);
			(*sample_f_data_y)(cur_row, sample_col_num - 1) = (*field_data[1])(2 * f_i + k, tmp_col);
			(*sample_index)(cur_row, sample_col_num - 1) = (*field_index)(2 * f_i + k, tmp_col);
			cur_row++;
		}
	}

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < sample_col_num - 1; j++) {
			int f_j = int(col*j / sample_col_num);
			int row_cnt = row - cur_row;
			(*sample_f_data_x)(cur_row, j) = (*field_data[0])(row - 2 + i, f_j);
			(*sample_f_data_y)(cur_row, j) = (*field_data[1])(row - 2 + i, f_j);
			(*sample_index)(cur_row, j) = (*field_index)(row - 2 + i, f_j);
		}
		if (cur_row % 2) tmp_col = col - 2;
		else tmp_col = col - 1;
		(*sample_f_data_x)(cur_row, sample_col_num - 1) = (*field_data[0])(row - 2 + i, tmp_col);
		(*sample_f_data_y)(cur_row, sample_col_num - 1) = (*field_data[1])(row - 2 + i, tmp_col);
		(*sample_index)(cur_row, sample_col_num - 1) = (*field_index)(row - 2 + i, tmp_col);
		cur_row++;
	}
	sample_field_data.push_back(sample_f_data_x);
	sample_field_data.push_back(sample_f_data_y);
	return sample_index;
}


MatrixXf* RectSdBkCalc::sample_field_data_pre(const int sample_row_num, const int sample_col_num)
{
	int row = solar_scene->layouts[0]->layout_row_col.x();
	int col = solar_scene->layouts[0]->layout_row_col.y();
	if (!sample_field_data.empty()) {
		delete sample_field_data[0];
		delete sample_field_data[1];
		sample_field_data.clear();
	}

	MatrixXf* sample_f_data_x = new MatrixXf(sample_row_num, sample_col_num);
	MatrixXf* sample_f_data_y = new MatrixXf(sample_row_num, sample_col_num);
	MatrixXf* sample_index = new MatrixXf(sample_row_num, sample_col_num);

	for (int i = 0; i < sample_row_num; i++) {
		int f_i = int(row*i / sample_row_num);
		for (int j = 0; j < sample_col_num; j++) {
			int f_j = int(col*j / sample_col_num);
			(*sample_f_data_x)(i, j) = (*field_data[0])(f_i, f_j);
			(*sample_f_data_y)(i, j) = (*field_data[1])(f_i, f_j);
			(*sample_index)(i, j) = (*field_index)(f_i, f_j);
		}
	}

	sample_field_data.push_back(sample_f_data_x);
	sample_field_data.push_back(sample_f_data_y);
	return sample_index;
}


void RectSdBkCalc::get_row_col(const int index, int& r, int & c)
{
	int col = solar_scene->layouts[0]->layout_row_col.y();
	r = index / col;
	c = index%col;
}


void FermatSdBkCalc::save_clipper_res(const string save_path, int month, int day, int hour, int minute)
{
	fstream outFile(save_path + "clipper_m" + to_string(month) + "_d" + to_string(day) + "_h" + to_string(hour) + "_min" + to_string(minute) + ".txt", ios_base::out);
	for (auto&h : solar_scene->helios) {
		outFile << h->helio_pos.x() << ' ' << h->helio_pos.z() << ' ' << h->sd_bk << endl;
	}
	outFile.close();
}


