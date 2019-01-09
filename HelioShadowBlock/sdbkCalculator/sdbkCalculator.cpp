#include"sdbkCalculator.h"
#include <random>
#include <ctime>


void SdBkCalc::sample_calc_preprocess(const int sample_row_num, const int sample_col_num, bool calc_s, bool calc_f)
{
	if (calc_f || field_index == nullptr)
		field_index = field_data_pre();
	if (calc_s || sample_field_index == nullptr)
		sample_field_index = sample_field_data_pre(sample_row_num, sample_col_num);
}

MatrixXd* SdBkCalc::calcSampleShadowBlock(const double DNI)
{
	vector<Receiver*> recvs = solar_scene->recvs;
	int row = sample_field_index->rows();
	int col = sample_field_index->cols();
	// Vector3d focus_center = recvs[0]->recv_pos + Vector3d(recvs[0]->recv_normal.array() * recvs[0]->recv_size.array());
	
	Vector3d reverse_sunray_dir = -solar_scene->sunray_dir;
	sample_sd_bk_res = new MatrixXd(row, col);

#pragma omp parallel for
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			int index = (*sample_field_index)(i, j);
			Heliostat* helio = solar_scene->helios[index];
			set<vector<int>> shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda;
			// Vector3d reflect_dir = focus_center - helio->helio_pos;
			int fc_index = solar_scene->helios[index]->focus_center_index;
			Vector3d reflect_dir = recvs[0]->focus_center[fc_index] - helio->helio_pos;
			calcIntersection3DDDA(helio, reverse_sunray_dir, shadow_relative_grid_label_3ddda);

			calcIntersection3DDDA(helio, reflect_dir, block_relative_grid_label_3ddda);
			vector<Vector3d> dir = { reverse_sunray_dir, reflect_dir };
			vector<set<vector<int>>> estimate_grids = { shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda };
			(*sample_sd_bk_res)(i, j) = helioClipper(helio, dir, estimate_grids);
			if (gl != NULL)
				helio->flux_sum = calcFluxMap(helio, DNI);
			//(*sample_sd_bk_res)(i, j) *= helio->flux_sum;

			shadow_relative_grid_label_3ddda.clear();
			block_relative_grid_label_3ddda.clear();
			estimate_grids.clear();
		}
	}
	return sample_sd_bk_res;
}


double SdBkCalc::helioClipper(Heliostat*helio, const Vector3d&dir, const set<vector<int>>& estimate_grid)
{
	vector<Vector3d> helio_v, local_v, tmp_v(4);
	vector<Vector2d> project_v(4);
	double t;
	Paths subj(1), clips;
	helio_v = helio->vertex;
	Vector3d reverse_dir = Vector3d(-dir.x(), -dir.y(), -dir.z());

	for (int i = 0; i < helio_v.size(); i++) {
		local_v.push_back(GeometryFunc::mulMatrix(helio_v[i], helio->world2localM));
		subj[0] << IntPoint(VERTEXSCALE*local_v[i].x(), VERTEXSCALE*local_v[i].z());
	}

	double total_area = 0;
	int subj_v_num = subj[0].size();
	for (int j = 0; j < subj_v_num; j++)
		total_area += (subj[0][j].X*subj[0][(j + 1) % subj_v_num].Y)
		- (subj[0][(j + 1) % subj_v_num].X*subj[0][j].Y);
	total_area = fabs(total_area*0.5);

	if (total_area == 0) {
		cout << "Project surface is 0!" << endl;
		return 0;
	}

	vector<Vector3d> pro(4);
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

	double sum = 0;
	for (int i = 0; i < solution.size(); i++) {
		int n = solution[i].size();
		for (int j = 0; j < n; j++) {
			sum += (solution[i][j].X*solution[i][(j + 1) % n].Y)
				- (solution[i][(j + 1) % n].X*solution[i][j].Y);
		}
	}
	sum = fabs(sum*0.5);

	double res = sum / total_area;

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

double SdBkCalc::helioClipper(Heliostat * helio, const vector<Vector3d>& dir, const vector<set<vector<int>>>& estimate_grids)
{
#ifdef CALC_TIME
	auto start = std::chrono::high_resolution_clock::now();

#endif // CALC_TIME


	vector<Vector3d> helio_v, tmp_v(4);
	vector<Vector2d> project_v(4), local_v;
	double t;
	Paths subj(1), clips;
	helio_v = helio->vertex;

	for (int i = 0; i < helio_v.size(); i++) {
		Vector3d tmp_v = GeometryFunc::mulMatrix(helio_v[i], helio->world2localM);
		local_v.push_back(Vector2d(tmp_v.x(), tmp_v.z()));
		subj[0] << IntPoint(VERTEXSCALE*local_v[i].x(), VERTEXSCALE*local_v[i].y());
	}

	double total_area = 0;
	// int subj_v_num = subj[0].size();
	for (int j = 0; j < local_v.size(); j++)
		total_area += (local_v[j].x()*local_v[(j + 1) % 4].y()
		- (local_v[(j + 1) % 4].x()*local_v[j].y()));
	total_area = fabs(total_area*0.5);

	if (total_area == 0) {
		cout << "Project surface is 0!" << endl;
		return 0;
	}

	helio->max_rela_dis = INT_MIN;
	helio->min_rela_dis = INT_MAX;

	for (int index = 0; index < 2; index++) {
		Vector3d reverse_dir = -dir[index]; 
		vector<Vector3d> pro(4), tmp_pro(4);
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
					//t = calcIntersectionPoint(helio_v[i], reverse_dir, helio->vertex[0], helio->vertex[1], helio->vertex[2]);
					t = GeometryFunc::calcIntersection(helio->helio_normal, helio->helio_pos, helio_v[i], reverse_dir, pro[i]);
					if(t > Epsilon){
						cnt++;
						if (helio->max_rela_dis < t)
							helio->max_rela_dis = t;
						if (helio->min_rela_dis > t)
							helio->min_rela_dis = t;
						//cout << relative_helio->helio_index << endl;
					}
				}
				if (cnt > 0) {
					//double rela_dis = (helio->helio_pos - relative_helio->helio_pos).norm();
					//if (helio->rela_dis < rela_dis)
					//	helio->rela_dis = rela_dis;

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

	double sum = 0;
	for (int i = 0; i < solution.size(); i++) {
		int n = solution[i].size();
		for (int j = 0; j < n; j++) {
			sum += (solution[i][j].X/(double)VERTEXSCALE*solution[i][(j + 1) % n].Y/(double)VERTEXSCALE)
				- (solution[i][(j + 1) % n].X/(double)VERTEXSCALE*solution[i][j].Y/(double)VERTEXSCALE);
		}
	}
	sum = fabs(sum*0.5);

	double res = sum / total_area;



#ifdef OUTPUTRES
	fstream outFile;
	outFile.open("sub_clipper.txt", ios_base::out);
	//outFile << "subj:" << endl;
	for (auto&tmp : subj[0])
		outFile << tmp.X / (double)VERTEXSCALE << ' ' << tmp.Y / (double)VERTEXSCALE << endl;
	outFile.close();
	//outFile << "clips:" << endl;
	outFile.open("clips.txt", ios_base::out);
	for (int i = 0; i < clips.size(); i++)
		for (int j = 0; j < clips[i].size(); j++)
			outFile << clips[i][j].X / (double)VERTEXSCALE << ' ' << clips[i][j].Y / (double)VERTEXSCALE << endl;
	outFile.close();
	//outFile << "solution:" << endl;
	outFile.open("solution.txt", ios_base::out);
	for (int i = 0; i < solution.size(); i++)
		for (int j = 0; j < solution[i].size(); j++)
			outFile << solution[i][j].X / (double)VERTEXSCALE << ' ' << solution[i][j].Y / (double)VERTEXSCALE << endl;
	outFile.close();

#endif // OUTPUTRES

#ifdef CALC_TIME
	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	std::cout << "clipper time: " << time << "s." << endl;

#endif // CALC_TIME


	return res;
}


//	使用ray tracing方式计算blocking与shadowing
//	将定日镜均匀剖分，每个小区域中心点 发出一根光线，与周围定日镜求交
//	记录产生blocking和shadowing的定日镜的编号
//	version: CPU
double SdBkCalc::calcAccurateIntersection(Heliostat * helio, const Vector3d & dir, set<vector<int>>& relative_helio_label)
{
	int helio_inner_rows = 20;
	int helio_inner_cols = 20;

	double total_sum = 0;
	int cnt = 0;

	vector<Vector3d> subhelio_v;
	helio->getSubHelioVertex(subhelio_v);
	int subHelio_num = subhelio_v.size() / 4;
	for (int i = 0; i < subHelio_num; i++) {
		vector<Vector3d> helio_v;
		for (int k = 0; k < 4; k++)
			helio_v.push_back(subhelio_v[4 * i + k]);
		Vector3d row_dir = (helio_v[1] - helio_v[0]) / helio_inner_rows;
		Vector3d col_dir = (helio_v[3] - helio_v[0]) / helio_inner_cols;

		for (int i = 0; i <= helio_inner_rows; i++) {
			for (int j = 0; j <= helio_inner_cols; j++) {
				total_sum++;
				Vector3d ori_v = helio_v[0] + i*row_dir + j*col_dir;
				double tMin = INT_MAX;
				Heliostat* hNear = nullptr;
				for (auto&h : solar_scene->helios) {
					if (h != helio) {
						vector<Vector3d> neigh_v = h->vertex;

						Vector3d E1 = neigh_v[1] - neigh_v[0];
						Vector3d E2 = neigh_v[2] - neigh_v[0];
						Vector3d pvec = dir.cross(E2);
						double det = E1.dot(pvec);

						// ray and triangle are parallel if det is close to 0
						if (fabsf(det) < Epsilon) continue;

						double invDet = 1 / det;

						Vector3d T = ori_v - neigh_v[0];
						Vector3d qvec = T.cross(E1);

						double t = E2.dot(qvec)*invDet;
						if (t < Epsilon) continue;

						Vector3d intersect_v = ori_v + dir*t;
						int l;
						for (l = 0; l < 4; l++) {
							Vector3d edg = neigh_v[(l + 1) % 4] - neigh_v[l];
							Vector3d line = intersect_v - neigh_v[l];
							Vector3d tmp_n = edg.cross(line);
							if (tmp_n.dot(Vector3d(0, 1, 0)) < Epsilon)
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
	double res = cnt / total_sum;

	return res;
}


double SdBkCalc::calcAccurateIntersection(Heliostat* helio, const vector<Vector3d>& dir, vector<set<vector<int>>>& relative_helio_label)
{
	int helio_inner_rows = 50;
	int helio_inner_cols = 50;

	double total_sum = 0;
	int cnt = 0;

	//vector<Vector3d> subhelio_v;
	//helio->getSubHelioVertex(subhelio_v);
	//int subHelio_num = subhelio_v.size() / 4;
	fstream outFile("ray" + to_string(helio->helio_index) + ".txt", ios_base::out);
	vector<set<Heliostat*>> sd_bk_h(2);
	for (int i = 0; i < helio->vertex.size(); i++) {
		vector<Vector3d> helio_v = helio->vertex;
		//vector<Vector3d> helio_v;
		//for (int k = 0; k < 4; k++)
		//	helio_v.push_back(subhelio_v[4 * i + k]);
		Vector3d row_dir = (helio_v[1] - helio_v[0]) / helio_inner_rows;
		Vector3d col_dir = (helio_v[3] - helio_v[0]) / helio_inner_cols;

		for (int i = 0; i <= helio_inner_rows; i++) {
			for (int j = 0; j <= helio_inner_cols; j++) {
				total_sum++;
				Vector3d ori_v = helio_v[0] + i*row_dir + j*col_dir;
				Vector3d proj_ori_v = GeometryFunc::mulMatrix(ori_v, helio->world2localM);
				outFile << proj_ori_v.x() << ' ' << proj_ori_v.z() << ' ';
				double tMin = INT_MAX;
				Heliostat* hNear = nullptr;
				bool sd = false;
				for (auto&h : solar_scene->helios) {
					if (h != helio) {
						if ((h->helio_pos - helio->helio_pos).norm() > 80)
							continue;

						vector<Vector3d> neigh_v = h->vertex;

						Vector3d E1 = neigh_v[1] - neigh_v[0];
						Vector3d E2 = neigh_v[2] - neigh_v[0];
						Vector3d pvec = dir[0].cross(E2);
						double det = E1.dot(pvec);

						// ray and triangle are parallel if det is close to 0
						if (fabsf(det) < Epsilon) continue;

						double invDet = 1 / det;

						Vector3d T = ori_v - neigh_v[0];
						Vector3d qvec = T.cross(E1);

						double t = E2.dot(qvec)*invDet;
						if (t < Epsilon) continue;

						Vector3d intersect_v = ori_v + dir[0]*t;
						int l;
						for (l = 0; l < 4; l++) {
							Vector3d edg = neigh_v[(l + 1) % 4] - neigh_v[l];
							Vector3d line = intersect_v - neigh_v[l];
							Vector3d tmp_n = edg.cross(line);
							if (tmp_n.dot(Vector3d(0, 1, 0)) < Epsilon)
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
					outFile << 1 << endl;
					sd_bk_h[0].insert(hNear);
					if (helio->max_rela_dis < tMin)
						helio->max_rela_dis = tMin;
				}
				else {
					for (auto&h : solar_scene->helios) {
						
						if (h != helio) {
							if ((h->helio_pos - helio->helio_pos).norm() > 80)
								continue;
		
							vector<Vector3d> neigh_v = h->vertex;

							Vector3d E1 = neigh_v[1] - neigh_v[0];
							Vector3d E2 = neigh_v[2] - neigh_v[0];
							Vector3d pvec = dir[1].cross(E2);
							double det = E1.dot(pvec);

							// ray and triangle are parallel if det is close to 0
							if (fabsf(det) < Epsilon) continue;

							double invDet = 1 / det;

							Vector3d T = ori_v - neigh_v[0];
							Vector3d qvec = T.cross(E1);

							double t = E2.dot(qvec)*invDet;
							if (t < Epsilon) continue;

							Vector3d intersect_v = ori_v + dir[1] * t;
							int l;
							for (l = 0; l < 4; l++) {
								Vector3d edg = neigh_v[(l + 1) % 4] - neigh_v[l];
								Vector3d line = intersect_v - neigh_v[l];
								Vector3d tmp_n = edg.cross(line);
								if (tmp_n.dot(Vector3d(0, 1, 0)) < Epsilon)
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
						outFile << 1 << endl;
						sd_bk_h[1].insert(hNear);
						if (helio->max_rela_dis < tMin)
							helio->max_rela_dis = tMin;
					}
					else
						outFile << 0 << endl;
				}
			}
		}
	}
	double res = cnt / total_sum;
	outFile.close();
	for (int i = 0; i < sd_bk_h.size(); i++) {
		cout << endl;
		for (auto& iter : sd_bk_h[i]) {
			cout << iter->helio_index << endl;
		}
	}
	return res;
}


// 使用3DDDA计算产生阴影和遮挡的相关定日镜
// version: CPU
void SdBkCalc::calcIntersection3DDDA(Heliostat * helio, const Vector3d & dir, set<vector<int>>& relative_helio_label)
{
#ifdef CALC_TIME
	auto start = std::chrono::high_resolution_clock::now();

#endif // CALC_TIME


	vector<Layout*>& layouts = solar_scene->layouts;
	vector<Vector3d> helio_v;
	vector<Vector2d> project_helio_v;

	helio_v = helio->vertex;
	for (int i = 0; i < 4; i++) {
		project_helio_v.push_back(Vector2d(helio_v[i].x(), helio_v[i].z()));
	}

	// 1. 確定光線在場地中y方向移動距離及最終離開網格的位置
	double upper_y = layouts[0]->layout_size.y() + layouts[0]->layout_bound_pos.y();
	Vector2d upper_v[4];
	for (int i = 0; i < 4; i++) {
		double dis = (upper_y - helio_v[i].y()) / dir.y();
		upper_v[i] = Vector2d(dis*dir.x(), dis*dir.z()) + project_helio_v[i];
	}

	// 2. 確定定日鏡反射光柱在鏡場中的矩形包圍盒範圍
	vector<Vector2d> boundBox(2);
	boundBox[0] = Vector2d(INT_MAX, INT_MAX);		// min boundary
	boundBox[1] = Vector2d(INT_MIN, INT_MIN);		// max boundary
	for (int i = 0; i < 4; i++) {
		boundBox[0].x() = fmin(boundBox[0].x(), fmin(project_helio_v[i].x(), upper_v[i].x()));
		boundBox[0].y() = fmin(boundBox[0].y(), fmin(project_helio_v[i].y(), upper_v[i].y()));
		boundBox[1].x() = fmax(boundBox[1].x(), fmax(project_helio_v[i].x(), upper_v[i].x()));
		boundBox[1].y() = fmax(boundBox[1].y(), fmax(project_helio_v[i].y(), upper_v[i].y()));
	}
	Vector2d ray_dir = Vector2d(dir.x(), dir.z()).normalized();
	Vector2d deltaT;
	Vector2d t;
	Vector2d origs[2] = {
		Vector2d(INT_MAX,INT_MAX),		// min vertex
		Vector2d(INT_MIN,INT_MIN)		// max vertex
	};

	
	// 3. 確定layout下各網格的矩陣間隔
	Vector2d cellDimension = Vector2d(layouts[0]->helio_interval.x(), layouts[0]->helio_interval.z());

	double colLength = boundBox[1].x() - boundBox[0].x();
	double rowLength = boundBox[1].y() - boundBox[0].y();

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
		Vector2d o_grid(
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

		double tmp = 0;
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
#ifdef CALC_TIME
	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	std::cout << "3ddda time: " << time << "s." << endl;

#endif // CALC_TIME


}


double SdBkCalc::calcIntersectionPoint(const Vector3d & orig, const Vector3d & dir, const Vector3d & A, const Vector3d & B, const Vector3d & C)
{
	Vector3d E1 = B - A;
	Vector3d E2 = C - A;
	Vector3d pvec = dir.cross(E2);
	double det = E1.dot(pvec);

	// ray and triangle are parallel if det is close to 0
	if (fabsf(det) < Epsilon) return 0;

	double invDet = 1 / det;

	Vector3d T = orig - A;

	Vector3d qvec = T.cross(E1);

	double t = E2.dot(qvec)*invDet;
	return t;
}

 
// helio: 待计算定日镜
// dir: 定日镜到接收器的反射光线方向
double SdBkCalc::calcFluxMap(Heliostat * helio, const double DNI)
{
#ifdef CALC_TIME
	auto start = std::chrono::high_resolution_clock::now();
#endif // CALC_TIME


	int fc_index = helio->focus_center_index;
	vector<Receiver*> recvs = solar_scene->recvs;
	Vector3d focus_center = recvs[0]->focus_center[fc_index];
	Vector3d reverse_dir = (helio->helio_pos - focus_center).normalized();		// The normal of image plane
	double _flux_sum1 = 0;
	double _flux_sum2 = 0;
	double _flux_sum3 = 0;
	double _flux_sum4 = 0;

	Matrix4d world2localM, local2worldM;
	GeometryFunc::getImgPlaneMatrixs(reverse_dir, focus_center, local2worldM, world2localM, 1);

	for (int i = 0; i < helio->cos_phi.size(); i++) {
		if (helio->cos_phi[i] > Epsilon) {
			vector<Vector2d> proj_v;
			vector<Vector3d> tmp_v;
			for (auto& v : recvs[0]->recv_vertex[i]) {
				Vector3d inter_v;
				GeometryFunc::calcIntersection(reverse_dir, focus_center, v, reverse_dir, inter_v);
				tmp_v.push_back(inter_v);
				inter_v = GeometryFunc::mulMatrix(inter_v, world2localM);
				proj_v.push_back(Vector2d(inter_v.x(), inter_v.z()));
			
			}
				//double _flux_sum_grid1 = _calc_flux_sum(proj_v, recvs[0]->mask_rows, recvs[0]->mask_cols, helio, helio->cos_phi[i], DNI);
				// double _flux_sum_grid2 = _calc_flux_sum(proj_v, helio, helio->cos_phi[i], DNI);
				double _flux_sum_grid3 = _multi_inte_flux_sum(proj_v, 2, helio, helio->cos_phi[i], DNI);
				// double _flux_sum_ray = ray_tracing_flux_sum(recvs[0]->recv_vertex[i], recvs[0]->focus_center[i], recvs[0]->recv_normal_list[i], helio, -reverse_dir, DNI);
				//_flux_sum1 += _flux_sum_grid1;
				// _flux_sum2 += _flux_sum_grid2;
				_flux_sum3 += _flux_sum_grid3;
				// _flux_sum4 += _flux_sum_ray;

				// flux_sum_matrix_grid(recvs[0]->recv_vertex[i], proj_v, recvs[0]->mask_rows, recvs[0]->mask_cols, helio, helio->cos_phi[i], DNI);
				// flux_sum_matrix_inte(recvs[0]->recv_normal_list[i], recvs[0]->focus_center[i], recvs[0]->recv_vertex[i], local2worldM, proj_v, helio, helio->cos_phi[i], DNI);
				// cout << _flux_sum1 << ' ' << _flux_sum2 << ' ' << _flux_sum3 << ' ' << _flux_sum4 << endl;
			
		}
	}
	//cout << _flux_sum1 << ' ' <<_flux_sum2 << ' ' << _flux_sum3 << endl;
#ifdef CALC_TIME
	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	std::cout << "Calculate flux sum time: " << time << "s." << endl;

#endif // CALC_TIME

	//fstream outFile("flux_sum.txt", ios_base::app);
	//outFile << _flux_sum1 << ' ' << _flux_sum2 << ' ' << _flux_sum3 << endl;
	return _flux_sum3;
}


void SdBkCalc::flux_sum_matrix_grid(vector<Vector3d>& _recv_v, vector<Vector2d>& proj_v, const int rows, const int cols, Heliostat* helio, const double cos_phi, const double DNI) {
	vector<Vector2d> recv_v;
	if (abs(_recv_v[0].x() - _recv_v[1].x())<Epsilon && abs(_recv_v[0].x() - _recv_v[2].x())<Epsilon)
		for (auto&v : _recv_v)
			recv_v.push_back(Vector2d(v.z(), v.y()));

	if (abs(_recv_v[0].z() - _recv_v[1].z()) < Epsilon && abs(_recv_v[0].z() - _recv_v[2].z())<Epsilon)
		for (auto&v : _recv_v)
			recv_v.push_back(Vector2d(v.x(), v.y()));


	MatrixXd recv_x(rows, cols), recv_y(rows, cols);
	MatrixXd mask_x(rows, cols), mask_y(rows, cols);	// 取中点密度作为该grid的密度
	Vector2d row_gap = (proj_v[2] - proj_v[1]) / cols;
	Vector2d col_gap = (proj_v[0] - proj_v[1]) / rows;
	Vector2d recv_row_gap = (recv_v[2] - recv_v[1]) / cols;
	Vector2d recv_col_gap = (recv_v[0] - recv_v[1]) / rows;
	Vector2d _v1 = proj_v[1] + 0.5*(row_gap + col_gap);
	Vector2d _v0 = proj_v[0] + 0.5*(row_gap - col_gap);
	Vector2d _v2 = proj_v[2] + 0.5*(col_gap - row_gap);
	Vector2d recv_v1 = recv_v[1] + 0.5*(recv_row_gap + recv_col_gap);
	Vector2d recv_v0 = recv_v[0] + 0.5*(recv_row_gap - recv_col_gap);
	Vector2d recv_v2 = recv_v[2] + 0.5*(recv_row_gap - recv_col_gap);
	for (int i = 0; i < rows; i++) {
		double start_v = _v1.x() + i * col_gap.x();
		double end_v = _v2.x() + i * col_gap.x();
		mask_x.row(i).setLinSpaced(cols, start_v, end_v);

		start_v = recv_v1.x() + i * recv_col_gap.x();
		end_v = recv_v2.x() + i * recv_col_gap.x();
		recv_x.row(i).setLinSpaced(cols, start_v, end_v);
	}
	for (int i = 0; i < cols; i++) {
		double start_v = _v1.y() + i * row_gap.y();
		double end_v = _v0.y() + i* row_gap.y();
		mask_y.col(i).setLinSpaced(rows, start_v, end_v);

		start_v = recv_v1.y() + i * recv_row_gap.y();
		end_v = recv_v0.y() + i* recv_row_gap.y();
		recv_y.col(i).setLinSpaced(rows, start_v, end_v);
	}
	double grid_area = 0.0;
	for (int i = 0; i < proj_v.size(); i++)
		grid_area += proj_v[i].x()*proj_v[(i + 1) % 4].y() - proj_v[i].y()*proj_v[(i + 1) % 4].x();
	grid_area = fabs(grid_area / 2.0 / rows / cols);

	double sum = 0.0;
	fstream outFile("grid_" + to_string(helio->helio_index) + ".txt", ios_base::out);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			// outFile << recv_x(i, j) << ' ' << recv_y(i, j) << ' ' 
			outFile << mask_x(i, j) << ' ' << mask_y(i, j) << ' '
				<< DNI*cos_phi*helio->flux_param * grid_area* gl->flux_func(mask_x(i, j), mask_y(i, j), helio->sigma, helio->l_w_ratio) << endl;
			sum += DNI*cos_phi*helio->flux_param * grid_area* gl->flux_func(mask_x(i, j), mask_y(i, j), helio->sigma, helio->l_w_ratio);
		}
	}
	outFile.close();
	cout << "grid: " << sum << endl;
}


void SdBkCalc::flux_sum_matrix_inte(Vector3d& recv_normal, Vector3d& fc, vector<Vector3d>& _recv_v, Matrix4d& local2world, vector<Vector2d>& proj_v, Heliostat * helio, const double cos_phi, const double DNI) {
	vector<VectorXd> weight = gl->getW();
	vector<VectorXd> node = gl->getX();
	Vector4d x, y;
	Vector4d recv_x, recv_y;
	x(0) = proj_v[0].x();
	y(0) = proj_v[0].y();
	x(1) = proj_v[3].x();
	y(1) = proj_v[3].y();
	x(2) = proj_v[2].x();
	y(2) = proj_v[2].y();
	x(3) = proj_v[1].x();
	y(3) = proj_v[1].y();

	vector<Vector2d> recv_v;
	if (abs(_recv_v[0].x() - _recv_v[1].x())<Epsilon && abs(_recv_v[0].x() - _recv_v[2].x())<Epsilon)
		for (auto&v : _recv_v)
			recv_v.push_back(Vector2d(v.z(), v.y()));

	if (abs(_recv_v[0].z() - _recv_v[1].z()) < Epsilon && abs(_recv_v[0].z() - _recv_v[2].z())<Epsilon)
		for (auto&v : _recv_v)
			recv_v.push_back(Vector2d(v.x(), v.y()));

	recv_x(0) = recv_v[0].x();
	recv_y(0) = recv_v[0].y();
	recv_x(1) = recv_v[3].x();
	recv_y(1) = recv_v[3].y();
	recv_x(2) = recv_v[2].x();
	recv_y(2) = recv_v[2].y();
	recv_x(3) = recv_v[1].x();
	recv_y(3) = recv_v[1].y();

	double sum = 0.0;
	fstream outFile("gauss_" + to_string(helio->helio_index) +".txt", ios_base::out);
	Vector2d map_v;
	Vector2d recv_map_v;
	for (int i = 0; i < weight[0].size(); i++) {
		for (int j = 0; j < weight[1].size(); j++) {
			map_v = gl->map(x, y, node[0][i], node[1][j]);
			recv_map_v = gl->map(recv_x, recv_y, node[0][i], node[1][j]);
			double tmp_sum = DNI*cos_phi*helio->flux_param *
				weight[0][i] * weight[1][j] * gl->jacobi(x, y, node[0](i), node[1](j))*gl->flux_func(map_v.x(), map_v.y(), helio->sigma, helio->l_w_ratio);
			Vector3d recv_map_v = GeometryFunc::mulMatrix(Vector3d(map_v.x(), 0, map_v.y()), local2world);
			// outFile << recv_map_v.x() << ' ' << recv_map_v.y() << ' ' << tmp_sum << endl;
			sum += tmp_sum;
			outFile << map_v.x() << ' ' << map_v.y() << ' ' << tmp_sum << endl;
		}
	}

	outFile.close();
	cout << "gauss: " << sum << endl;
}

///
// 对接收器投影面进行均匀分割，统计每个分割面的flux分布结果，并求和
///
double SdBkCalc::_calc_flux_sum(vector<Vector2d>& proj_v, const int rows, const int cols, Heliostat * helio, const double cos_phi, const double DNI)
{
	auto start = std::chrono::high_resolution_clock::now();

	MatrixXd mask_x(rows, cols), mask_y(rows, cols), tmp_mask_y(rows, cols);	// 取中点密度作为该grid的密度
	Vector2d row_gap = (proj_v[2] - proj_v[1]) / cols;
	Vector2d col_gap = (proj_v[0] - proj_v[1]) / rows;
	Vector2d _v1 = proj_v[1] + 0.5*(row_gap + col_gap);
	Vector2d _v0 = proj_v[0] + 0.5*(row_gap - col_gap);
	Vector2d _v2 = proj_v[2] + 0.5*(col_gap - row_gap);
	for (int i = 0; i < rows; i++) {
		double start_v = _v1.x() + i * col_gap.x();
		double end_v = _v2.x() + i * col_gap.x();
		mask_x.row(i).setLinSpaced(cols, start_v, end_v);
	}
	for (int i = 0; i < cols; i++) {
		double start_v = _v1.y() + i * row_gap.y();
		double end_v = _v0.y() + i* row_gap.y();
		tmp_mask_y.col(i).setLinSpaced(rows, start_v, end_v);
		mask_y.col(i).setLinSpaced(rows, helio->l_w_ratio * start_v, helio->l_w_ratio * end_v);
	}

	double flux_sum = (-0.5 / pow(helio->sigma, 2) * (mask_x.array().pow(2) + mask_y.array().pow(2))).array().exp().sum();
	double grid_area = 0.0;
	for (int i = 0; i < proj_v.size(); i++)
		grid_area += proj_v[i].x()*proj_v[(i + 1) % 4].y() - proj_v[i].y()*proj_v[(i + 1) % 4].x();
	grid_area = fabs(grid_area / 2.0 / rows / cols);
	// flux_sum = flux_sum * (1 - helio->sd_bk)*DNI*cos_phi*helio->flux_param * grid_area;
	flux_sum = flux_sum *DNI*cos_phi*helio->flux_param * grid_area;
	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	// std::cout << "Calculate flux sum time: " << time << "s." << endl;

	return flux_sum;
}


double SdBkCalc::_calc_flux_sum(vector<Vector2d>& proj_v, Heliostat * helio, const double cos_phi, const double DNI)
{
	vector<VectorXd> weight = gl->getW();
	vector<VectorXd> node = gl->getX();
	Vector4d x, y;
	x(0) = proj_v[0].x();
	y(0) = proj_v[0].y();
	x(1) = proj_v[3].x();
	y(1) = proj_v[3].y();
	x(2) = proj_v[2].x();
	y(2) = proj_v[2].y();
	x(3) = proj_v[1].x();
	y(3) = proj_v[1].y();

	double sum = gl->calcInte(x, y, helio->sigma, helio->l_w_ratio);
	//sum = sum * helio->flux_param * (1 - helio->sd_bk) * DNI * cos_phi;
	sum = sum * helio->flux_param *  DNI * cos_phi;

	return sum;
}


double SdBkCalc::_multi_inte_flux_sum(vector<Vector2d>& proj_v, int n, Heliostat* helio, const double cos_phi, const double DNI) {
	Vector2d row_gap = (proj_v[3] - proj_v[0]) / n;
	Vector2d col_gap = (proj_v[1] - proj_v[0]) / n;

	Vector4d tmp_x, tmp_y;
	double sum = 0.0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			tmp_x(0) = (proj_v[0] + i*row_gap + j*col_gap).x();
			tmp_y(0) = (proj_v[0] + i*row_gap + j*col_gap).y();
			tmp_x(1) = (proj_v[0] + (i + 1)*row_gap + j*col_gap).x();
			tmp_y(1) = (proj_v[0] + (i + 1)*row_gap + j*col_gap).y();
			tmp_x(2) = (proj_v[0] + (i + 1)*row_gap + (j + 1)*col_gap).x();
			tmp_y(2) = (proj_v[0] + (i + 1)*row_gap + (j + 1)*col_gap).y();
			tmp_x(3) = (proj_v[0] + i*row_gap + (j + 1)*col_gap).x();
			tmp_y(3) = (proj_v[0] + i*row_gap + (j + 1)*col_gap).y();
			sum += gl->calcInte(tmp_x, tmp_y, helio->sigma, helio->l_w_ratio);
		}
	}	
	//sum = sum * helio->flux_param * (1 - helio->sd_bk) * DNI * cos_phi;
	sum = sum * helio->flux_param *  DNI * cos_phi;
	return sum;
}


inline double calc_mAA(double dis) {
	double mAA;
	if (dis <= 1000)
		mAA = (double)(0.99321 - 0.0001176 * dis + 1.97 * 1e-8 * dis * dis);      //d<1000
	else
		mAA = exp(-0.0001106 * dis);
	return mAA;
}


inline double random_double() {
	static std::default_random_engine e(time(NULL));
	static std::uniform_real_distribution<double> u(0, 1);
	return u(e);
}

///
// 使用Ray-Tracing计算接收器上接收到的能量
///
double SdBkCalc::ray_tracing_flux_sum(vector<Vector3d>& recv_v, Vector3d& recv_pos, Vector3d& recv_normal, Heliostat * helio, const Vector3d& dir, const double DNI)
{
	 int rows = helio->helio_size.z() / RECEIVER_SLICE;
	 int cols = helio->helio_size.x() / RECEIVER_SLICE;
	 int grid_num = rows*cols;
	int ray_num = 20;
	double ray_e = DNI*helio->S * helio->cos_w / grid_num / ray_num;
	int iter_n = 2;
	double sum =0;
	Vector3d row_gap = (helio->vertex[3] - helio->vertex[0]) / cols;
	Vector3d col_gap = (helio->vertex[1] - helio->vertex[0]) / rows;

	Vector3d start_v;
	Vector3d ray_v, inter_v, edg, line, tmp_n;
	for (int n = 0; n < iter_n; n++) {
		start_v = helio->vertex[0];
		for (int i = 0; i < rows; i++) {
			start_v = helio->vertex[0] + i*col_gap;
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < ray_num; k++) {
					double x = random_double();
					double y = random_double();
					ray_v = start_v + x*row_gap + y*col_gap;
					if (GeometryFunc::calcIntersection(recv_normal, recv_pos, ray_v, dir, inter_v) < Epsilon && GeometryFunc::inProjArea(recv_v, inter_v)) {
						//int l = 0;
						//for (; l < 4; l++) {
						//	edg = recv_v[(l + 1) % 4] - recv_v[l];
						//	line = inter_v - recv_v[l];
						//	tmp_n = edg.cross(line);
						//	if (tmp_n.dot(recv_normal) < Epsilon)
						//		break;
						//}
						//if (l == 4) {
						double dis = (inter_v - ray_v).norm();
						double mAA = calc_mAA(dis);
						sum += mAA*ray_e*helio->rou;
						//}
					}
				}
				start_v += row_gap;
			}
		}
	}

	return sum / iter_n;
}


///
// 数值积分在无穷大平面上积分得到结果
///
double SdBkCalc::inte_infinite_flux_sum(Heliostat * helio, const Vector3d& recv_pos,  const double cos_phi, const double DNI)
{
	double dis = (helio->helio_pos - recv_pos).norm();
	double mAA = calc_mAA(dis);
	return DNI *helio->S *helio->cos_w * helio->rou* cos_phi* mAA;
}


// 以光线跟踪计算结果为groundtruth，检测预测结果是否包含准确结果
// version: CPU
double SdBkCalc::checkForRelativeHelio(const set<vector<int>>& accurate_helio, const set<vector<int>>& estimate_helio)
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
	//return sum > 0 ? (double)cnt / (double)sum : 0;
	return (double)cnt / (double)solar_scene->helios.size();
}


double SdBkCalc::calcSingleShadowBlock(int helio_index)
{
	vector<Receiver*> recvs = solar_scene->recvs;

	auto helio = this->solar_scene->helios[helio_index];
	set<vector<int>> shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda;
	Vector3d reverse_sunray_dir = -solar_scene->sunray_dir;
	// Vector3d focus_center = recvs[0]->recv_pos + Vector3d(recvs[0]->recv_normal.array() * recvs[0]->recv_size.array());
	int fc_index = solar_scene->helios[helio_index]->focus_center_index;
	Vector3d reflect_dir = recvs[0]->focus_center[fc_index] - helio->helio_pos;

	// 3DDDA + clipper 
	//Calc the relative heliostats which cause shadowing
	calcIntersection3DDDA(helio, reverse_sunray_dir, shadow_relative_grid_label_3ddda);

	//Calc the relactive heliostats which cause blocking
	// Vector3d reflect_dir = (focus_center - helio->helio_pos).normalized();
	calcIntersection3DDDA(helio, reflect_dir, block_relative_grid_label_3ddda);

	vector<Vector3d> dir = { reverse_sunray_dir, reflect_dir };
	vector<set<vector<int>>> estimate_grids = { shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda };
	return helioClipper(helio, dir, estimate_grids);
}


double SdBkCalc::calcSingleFluxSum(int helio_index, const double DNI) {

	Heliostat* helio = solar_scene->helios[helio_index];
	int fc_index = helio->focus_center_index;
	vector<Receiver*> recvs = solar_scene->recvs;
	Vector3d focus_center = recvs[0]->focus_center[fc_index];
	Vector3d reverse_dir = (helio->helio_pos - focus_center).normalized();		// The normal of image plane
	double _flux_sum1 = 0;
	double _flux_sum2 = 0;
	double _flux_sum3 = 0;
	double _flux_sum4 = 0;

	Matrix4d world2localM, local2worldM;
	GeometryFunc::getImgPlaneMatrixs(reverse_dir, focus_center, local2worldM, world2localM, 1);

	for (int i = 0; i < helio->cos_phi.size(); i++) {
		if (helio->cos_phi[i] > Epsilon) {
			vector<Vector2d> proj_v;
			vector<Vector3d> tmp_v;
			for (auto& v : recvs[0]->recv_vertex[i]) {
				Vector3d inter_v;
				GeometryFunc::calcIntersection(reverse_dir, focus_center, v, reverse_dir, inter_v);
				tmp_v.push_back(inter_v);
				inter_v = GeometryFunc::mulMatrix(inter_v, world2localM);
				proj_v.push_back(Vector2d(inter_v.x(), inter_v.z()));

			}
			//double _flux_sum_grid1 = _calc_flux_sum(proj_v, recvs[0]->mask_rows, recvs[0]->mask_cols, helio, helio->cos_phi[i], DNI);
			// double _flux_sum_grid2 = _calc_flux_sum(proj_v, helio, helio->cos_phi[i], DNI);
			double _flux_sum_grid3 = _multi_inte_flux_sum(proj_v, 2, helio, helio->cos_phi[i], DNI);
			// double _flux_sum_ray = ray_tracing_flux_sum(recvs[0]->recv_vertex[i], recvs[0]->focus_center[i], recvs[0]->recv_normal_list[i], helio, -reverse_dir, DNI);
			//_flux_sum1 += _flux_sum_grid1;
			// _flux_sum2 += _flux_sum_grid2;
			_flux_sum3 += _flux_sum_grid3;
			// _flux_sum4 += _flux_sum_ray;

		}
	}

	return _flux_sum3;
}


void SdBkCalc::calcShadowBlock(const double DNI)
{
	vector<Heliostat*> helios = solar_scene->helios;
	vector<Layout*> layouts = solar_scene->layouts;
	vector<Receiver*> recvs = solar_scene->recvs;

	//int row = layouts[0]->layout_row_col.x();
	//int col = layouts[0]->layout_row_col.y();
	//sd_bk_res = new MatrixXd(row, col);

//#ifdef DEBUG
//	MatrixXd* raytracing_store = new MatrixXd(row, col);
//#endif // DEBUG


#ifdef READFILE
	fstream inFile("shadowblock_gt_save.txt", ios_base::in);
	if (inFile.fail())
		cout << "Can't open the file!" << endl;
#endif // READFILE

#ifdef CALC_TIME
	auto start = std::chrono::high_resolution_clock::now();

#endif // CALC_TIME

	Vector3d reverse_sunray_dir = -solar_scene->sunray_dir;
#pragma omp parallel for
	for (int i = 0; i < helios.size(); i++) {

		auto helio = helios[i];
		set<vector<int>> shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda;

		// 3DDDA + clipper 
		//Calc the relative heliostats which cause shadowing
		calcIntersection3DDDA(helio, reverse_sunray_dir, shadow_relative_grid_label_3ddda);

		//Calc the relactive heliostats which cause blocking
		int fc_index = helio->focus_center_index;
		Vector3d reflect_dir = (recvs[0]->focus_center[fc_index] - helio->helio_pos).normalized();

		calcIntersection3DDDA(helio, reflect_dir, block_relative_grid_label_3ddda);

		vector<Vector3d> dir = { reverse_sunray_dir, reflect_dir };		// from heliostat
		if (shadow_relative_grid_label_3ddda.size() != 0 || block_relative_grid_label_3ddda.size() != 0) {
			vector<set<vector<int>>> estimate_grids = { shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda };
			helio->sd_bk = helioClipper(helio, dir, estimate_grids);
		}
		else
			helio->sd_bk = 0;

		if(gl!=NULL)
			helio->flux_sum = calcFluxMap(helio, DNI);

		double tanpi2zen = reverse_sunray_dir.y() / sqrt(pow(reverse_sunray_dir.x(),2) + pow(reverse_sunray_dir.z(), 2));
		double HIh = max(helio->helio_size.x(), helio->helio_size.z());
		helio->approx_rela_dis = (HIh*sin(acos(helio->helio_normal.y()))) / tanpi2zen + HIh*helio->helio_normal.y();
		//cout << helio->sd_bk << endl;
#ifdef DEBUG
		// Ray tracing
		vector<set<vector<int>>> ray_relative_grids(2);
		int helio_col = (helio->helio_pos.x() - layouts[0]->layout_first_helio_center.x()) / layouts[0]->helio_interval.x();			// smaller x is, smaller col is
		int helio_row = (helio->helio_pos.z() - layouts[0]->layout_first_helio_center.z()) / layouts[0]->helio_interval.z();


		double ray_res = calcAccurateIntersection(helio, dir, ray_relative_grids);
		cout << "Heliostat: " << helio->helio_index << ' ' << ray_res << ' ' << helio->sd_bk << endl;		helio->sd_bk = ray_res;
		
		//int r1 = (helio->helio_pos.z() - solar_scene->layouts[0]->layout_first_helio_center.z()) / solar_scene->layouts[0]->helio_interval.z();		// smaller x is, smaller col is
		//int c1 = (helio->helio_pos.x() - solar_scene->layouts[0]->layout_first_helio_center.x()) / solar_scene->layouts[0]->helio_interval.x();		// bigger z is, smaller row is
		//cout << r1 << ' ' << c1 << endl;
		//(*raytracing_store)(r,c) = calcAccurateIntersection(helio, dir);
		//cout << (*raytracing_store)(r, c) - (*sd_bk_res)(r, c) << endl;
		for (auto&label : ray_relative_grids[0]) {
			if(shadow_relative_grid_label_3ddda.count(label) == 0)
				cout << helio_col << ' ' << helio_row << ' ' << helio->helio_index << endl;
		}
		for (auto&label : ray_relative_grids[1]) {
			if (block_relative_grid_label_3ddda.count(label) == 0) {
				cout << helio_col << ' ' << helio_row << ' ' << helio->helio_index << endl;
			}
		}
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

#ifdef CALC_TIME
	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	std::cout << "Heliostat total calculate time: " << time << "s." << endl;

#endif // CALC_TIME



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


MatrixXd* CrossRectSdBkCalc::field_data_pre() {
	vector<Heliostat*>& helios = solar_scene->helios;
	vector<Layout*>& layouts = solar_scene->layouts;
	int row = layouts[0]->layout_row_col.x();
	int col = layouts[0]->layout_row_col.y();
	if (!field_data.empty()) {
		delete field_data[0];
		delete field_data[1];
		field_data.clear();
	}

	MatrixXd* field_data_x = new MatrixXd(row, col);
	MatrixXd* field_data_y = new MatrixXd(row, col);
	MatrixXd* field_index = new MatrixXd(row, col);
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


MatrixXd* CrossRectSdBkCalc::sample_field_data_pre(const int sample_row_num, const int sample_col_num)
{
	int row = solar_scene->layouts[0]->layout_row_col.x();
	int col = solar_scene->layouts[0]->layout_row_col.y();
	if (!sample_field_data.empty()) {
		delete sample_field_data[0];
		delete sample_field_data[1];
		sample_field_data.clear();
	}

	MatrixXd* sample_f_data_x = new MatrixXd(sample_row_num, sample_col_num);
	MatrixXd* sample_f_data_y = new MatrixXd(sample_row_num, sample_col_num);
	MatrixXd* sample_index = new MatrixXd(sample_row_num, sample_col_num);
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


MatrixXd* RectSdBkCalc::sample_field_data_pre(const int sample_row_num, const int sample_col_num)
{
	int row = solar_scene->layouts[0]->layout_row_col.x();
	int col = solar_scene->layouts[0]->layout_row_col.y();
	if (!sample_field_data.empty()) {
		delete sample_field_data[0];
		delete sample_field_data[1];
		sample_field_data.clear();
	}

	MatrixXd* sample_f_data_x = new MatrixXd(sample_row_num, sample_col_num);
	MatrixXd* sample_f_data_y = new MatrixXd(sample_row_num, sample_col_num);
	MatrixXd* sample_index = new MatrixXd(sample_row_num, sample_col_num);

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
	for (int i = 0; i < solar_scene->helios.size(); i++) {
		auto h = solar_scene->helios[i];
		outFile << solar_scene->sunray_dir.x() << ' ' << solar_scene->sunray_dir.z() << ' ' << h->helio_pos.x() << ' ' << h->helio_pos.z() << ' ' << h->sd_bk << ' ' 
			<< h->flux_sum << ' ' << h->min_rela_dis << ' ' << h->max_rela_dis << ' ' << h->approx_rela_dis << endl;
	}
	outFile.close();

}


