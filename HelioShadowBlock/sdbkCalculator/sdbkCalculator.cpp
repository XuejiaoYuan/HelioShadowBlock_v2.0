#include"sdbkCalculator.h"
#include <random>
#include <ctime>
#include "../DataStructure/Timer.h"

//
// [计算采样定日镜能量] 计算每个采样定日镜反射到接收器上课被接收器吸收的能量
//
void SdBkCalc::calcSampleEnergy(int sample_row, int sample_col, const double DNI) {
	int helio_sum = solar_scene->helios.size();
	double gap = (double)helio_sum / (sample_row*sample_col);	

#pragma omp parallel for
	for(int i=0; i<sample_row; ++i)
		for (int j = 0; j < sample_col; ++j) {
			int index = (i*sample_col*gap + j*gap);
			_helio_calc(index, DNI);
		}
}

void SdBkCalc::saveCalcRes(const string s)
{
	fstream outFile(s, ios_base::out);
	for (auto&h : solar_scene->helios) {
		outFile << h->helio_pos.x() << ' ' << h->helio_pos.z() << ' ' << h->total_e << endl;
	}
	outFile.close();
}


//
// [计算单个定日镜能量] 返回单个定日镜反射到接收器上被接收器吸收的能量
//
double SdBkCalc::_helio_calc(int index, int DNI)
{
	Heliostat* helio = solar_scene->helios[index];
	set<vector<int>> shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda;
	int fc_index = solar_scene->helios[index]->focus_center_index;
	Vector3d reflect_dir = solar_scene->recvs[0]->focus_center[fc_index] - helio->helio_pos;
	calcIntersection3DDDA(helio, -solar_scene->sunray_dir, shadow_relative_grid_label_3ddda);

	calcIntersection3DDDA(helio, reflect_dir, block_relative_grid_label_3ddda);
	vector<Vector3d> dir = { -solar_scene->sunray_dir, reflect_dir };
	vector<set<vector<int>>> estimate_grids = { shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda };
	helio->sd_bk = helioClipper(helio, dir, estimate_grids);
	if (gl != NULL)
		helio->flux_sum = calcFluxMap(helio, DNI);
	helio->total_e = (1 - helio->sd_bk)*helio->flux_sum;
	helio->fluxCalc = true;

	shadow_relative_grid_label_3ddda.clear();
	block_relative_grid_label_3ddda.clear();
	estimate_grids.clear();
	return helio->total_e;
}


//
// [多边形裁剪] 单独处理阴影或遮挡相关定日镜
//
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

//
// [多边形裁剪] 处理阴影和遮挡，考虑阴影遮挡重叠部分
//
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


//
// [3DDDA] 
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


//
// [计算交点]
//
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


//
// [计算通量密度]
// helio: 待计算定日镜
// dir: 定日镜到接收器的反射光线方向
double SdBkCalc::calcFluxMap(Heliostat * helio, const double DNI)
{
#ifdef CALC_TIME
	Timer::resetStart();
#endif // CALC_TIME

	int fc_index = helio->focus_center_index;
	vector<Receiver*> recvs = solar_scene->recvs;
	Vector3d focus_center = recvs[0]->focus_center[fc_index];
	Vector3d reverse_dir = (helio->helio_pos - focus_center).normalized();		// The normal of image plane
	double _flux_sum = 0;

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
			_flux_sum += _multi_inte_flux_sum(proj_v, 2, helio, helio->cos_phi[i], DNI);;
		}
	}

#ifdef CALC_TIME
	Timer::printDuration("Calculate flux sum time");
#endif // CALC_TIME

	return _flux_sum;
}

//
// [采样计算通量密度] 以区域中心点通量密度代表该区域通量平均通量密度
//		计算每个点的结果并存入文件中
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
			outFile << mask_x(i, j) << ' ' << mask_y(i, j) << ' '
				<< DNI*cos_phi*helio->flux_param * grid_area* gl->flux_func(mask_x(i, j), mask_y(i, j), helio->sigma, helio->l_w_ratio) << endl;
			sum += DNI*cos_phi*helio->flux_param * grid_area* gl->flux_func(mask_x(i, j), mask_y(i, j), helio->sigma, helio->l_w_ratio);
		}
	}
	outFile.close();
	cout << "grid: " << sum << endl;
}


//
// [卷积计算通量密度] 以高斯积分的方式计算区域内的通量密度
//		计算每个点的结果并存入文件中
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
// [采样计算通量密度] 以区域中心点通量密度代表该区域通量平均通量密度
// 对接收器投影面进行均匀分割，统计每个分割面的flux分布结果，并求和
///
double SdBkCalc::_calc_flux_sum(vector<Vector2d>& proj_v, const int rows, const int cols, Heliostat * helio, const double cos_phi, const double DNI)
{
	Timer::resetStart();

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
	flux_sum = flux_sum *DNI*cos_phi*helio->flux_param * grid_area;
	
	Timer::printDuration("Calculate flux sum time");

	return flux_sum;
}


//
// [卷积计算通量密度] 以整体面积作为卷积
//
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
	sum *= helio->flux_param *  DNI * cos_phi;

	return sum;
}

//
// [卷积计算通量密度] 将区域分割成若干子区域，对每个子区域进行卷积
//
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
// [Ray Tracing计算能量] 使用Ray-Tracing计算接收器上接收到的能量
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

// 
// [检查3DDDA相关性是否可靠] 以光线跟踪计算结果为groundtruth，检测预测结果是否包含准确结果
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


//
// [计算定日镜的阴影遮挡] 计算目标定日镜的阴影与遮挡结果
//
double SdBkCalc::calcSingleShadowBlock(int helio_index)
{
	Heliostat* helio = solar_scene->helios[helio_index];
	set<vector<int>> shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda;
	int fc_index = solar_scene->helios[helio_index]->focus_center_index;
	Vector3d reflect_dir = solar_scene->recvs[0]->focus_center[fc_index] - helio->helio_pos;
	calcIntersection3DDDA(helio, -solar_scene->sunray_dir, shadow_relative_grid_label_3ddda);

	calcIntersection3DDDA(helio, reflect_dir, block_relative_grid_label_3ddda);
	vector<Vector3d> dir = { -solar_scene->sunray_dir, reflect_dir };
	vector<set<vector<int>>> estimate_grids = { shadow_relative_grid_label_3ddda, block_relative_grid_label_3ddda };
	helio->sd_bk = helioClipper(helio, dir, estimate_grids);

	return helio->sd_bk;
}

//
// [计算定日镜的通量密度和] 计算目标定日镜投射到接收器面上的能量和
//
double SdBkCalc::calcSingleFluxSum(int helio_index, const double DNI) {

	Heliostat* helio = solar_scene->helios[helio_index];
	int fc_index = helio->focus_center_index;
	vector<Receiver*> recvs = solar_scene->recvs;
	Vector3d focus_center = recvs[0]->focus_center[fc_index];
	Vector3d reverse_dir = (helio->helio_pos - focus_center).normalized();		// The normal of image plane
	double _flux_sum = 0;

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
			_flux_sum += _multi_inte_flux_sum(proj_v, 2, helio, helio->cos_phi[i], DNI);
		}
	}

	return _flux_sum;
}


//
// [计算所有定日镜反射能量] 计算所有定日镜反射到接收器上，接收器可获得的能量总和
//
double SdBkCalc::calcTotalEnergy(const double DNI)
{
	vector<Heliostat*> helios = solar_scene->helios;
	double sum = 0.0;
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
		double res  = _helio_calc(i, DNI);

#pragma omp critical
		sum += res;

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
	return  sum;
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


//
// [测试镜场所有定日镜阴影遮挡预处理函数]
//
void SdBkCalcTest::totalHeliosTest() {
	for (int i = 0; i < solar_scene->helios.size(); ++i) {
		singleHelioTest(i);
	}
}

//
// [测试单个定日镜阴影遮挡预处理函数]
//
void SdBkCalcTest::singleHelioTest(const int _helio_index) {
	setDir(-solar_scene->sunray_dir, solar_scene->helios[_helio_index]->helio_normal);
	setTestIndex(_helio_index);

	rayTracingSdBk();
	normalSdBk();
	boundingSphereSdBk();
	neighRowSdBk();
	improvedNeighSdBk();
	use3dddaSdBk();
}


// [ray tracing] 处理阴影或遮挡
//	使用ray tracing方式计算blocking和shadowing
//	将定日镜均匀剖分，每个小区域中心点 发出一根光线，与周围定日镜求交
//	记录产生blocking和shadowing的定日镜的编号
//	version: CPU
void SdBkCalcTest::rayTracingSdBk() {
	Timer::resetStart();

	int helio_inner_rows = 20;
	int helio_inner_cols = 20;

	double total_sum = helio_inner_rows *helio_inner_cols;
	int cnt = 0;

	auto helio = solar_scene->helios[helio_index];
	auto helio_v = helio->vertex;
	Vector3d row_dir = (helio_v[1] - helio_v[0]) / helio_inner_rows;
	Vector3d col_dir = (helio_v[3] - helio_v[0]) / helio_inner_cols;

#pragma omp parallel for
	for (int i = 0; i <= helio_inner_rows; i++) {
		for (int j = 0; j <= helio_inner_cols; j++) {
			Vector3d ori_v = helio_v[0] + i*row_dir + j*col_dir;
			if (calcIntersect(ori_v, sd_dir, gt_sd_helio_index)) {
			#pragma omp critical
				++cnt;
			}
			else if (calcIntersect(ori_v, bk_dir, gt_bk_helio_index)) {
			#pragma omp critical
				++cnt;
			}
				
		}
	}

	double res = cnt / total_sum;

	double time = Timer::getDuration();
	fstream outFile(save_path + "/rayTracing_time.txt", ios_base::app || ios_base::out);
	outFile << time << endl;
	outFile.close();
	outFile.open(save_path + "/rayTracing_sdbk.txt", ios_base::app || ios_base::out);
	outFile << res << endl;
	outFile.close();
	outFile.open(save_path + "/rayTracing_sd_index.txt", ios_base::app || ios_base::out);
	for (auto&n : gt_sd_helio_index)
		outFile << n << endl;
	outFile.close();
	outFile.open(save_path + "/rayTracing_bk_index.txt", ios_base::app || ios_base::out);
	for (auto&n : gt_bk_helio_index)
		outFile << n << endl;
	outFile.close();
}

//
// [ray tracing] 计算阴影或遮挡
//		返回当前起始点发射的光线是否与定日镜相交
bool SdBkCalcTest::calcIntersect(Vector3d& ori_v, Vector3d& dir, set<int>& index_set) {
	double tMin = INT_MAX;
	Heliostat* hNear = nullptr;

	for (int i = 0; i < solar_scene->helios.size(); ++i) {
		auto h = solar_scene->helios[i];
		if (h->helio_index != helio_index) {
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
		index_set.insert(hNear->helio_index);
		return true;
	}
	return false;
}


//
// [法向剔除] 计算两个定日镜之间的法向点积并进行剔除
//		设置相关定日镜的最远距离，以减少相关定日镜数量
void SdBkCalcTest::normalSdBk() {
	Heliostat* cur_h = solar_scene->helios[helio_index];
	set<int> sd_set, bk_set;
	int sd_ac = 0;
	int bk_ac = 0;
#pragma omp parallel for
	for (int i = 0; i < solar_scene->helios.size(); ++i) {
		auto h = solar_scene->helios[i];
		if (h->helio_index != helio_index) {
			if (h->helio_normal.dot(sd_dir) > Epsilon) {
				if (checkHelioDis(sd_dir, cur_h, cur_h->helio_pos, h->helio_pos)) {
					#pragma omp critical
					sd_set.insert(h->helio_index);
					if (gt_sd_helio_index.count(h->helio_index)) {
					#pragma omp critical
						++sd_ac;
					}
				}
			}
			if (h->helio_normal.dot(bk_dir) > Epsilon) {
				if (checkHelioDis(bk_dir, cur_h, cur_h->helio_pos, h->helio_pos)) {
					#pragma omp critical
					bk_set.insert(h->helio_index);
					if (gt_bk_helio_index.count(h->helio_index)) {
					#pragma omp critical
						++bk_ac;
					}
				}
			}
		}
	}

	saveTestRes("normal", sd_ac, bk_ac, sd_set, bk_set);
}


//
// [法向剔除] 法向剔除輔助函數，根據距離排除無關定日鏡
//
bool SdBkCalcTest::checkHelioDis(Vector3d& dir, Heliostat* helio, Vector3d& Hloc, Vector3d& HIloc) {
	double tanpi2zen = dir.y() / sqrt(dir.x()*dir.x() + dir.z()*dir.z());
	Vector3d HIt = helio->helio_normal;
	double HIzen = acos(HIt.y());
	double HIh = helio->helio_size.x();
	double l_max = (HIloc.y() - Hloc.y() + HIh*sin(HIzen)) / tanpi2zen + HIh*HIt.y();
	double interaction_limit = 100;
	l_max = fmin(l_max, interaction_limit*HIh);	//limit to a reasonable number
	double hdist = (HIloc - Hloc).norm();
	if (hdist < l_max) return true;
	else return false;
}

//
// [bounding sphere] 通過計算兩個定日鏡包圍盒之間的距离判断是否相关
//
void SdBkCalcTest::boundingSphereSdBk() {
	Heliostat* cur_h = solar_scene->helios[helio_index];
	set<int> sd_set, bk_set;
	int sd_ac = 0;
	int bk_ac = 0;
	double diameter = sqrt(pow(cur_h->helio_size.x(), 2) + pow(cur_h->helio_size.z(), 2));
#pragma omp parallel for
	for (int i = 0; i < solar_scene->helios.size(); ++i) {
		auto h = solar_scene->helios[i];
		if (h->helio_index != helio_index) {
			if(checkBoundingBox(cur_h->helio_pos, h->helio_pos, sd_dir, diameter)) {
#pragma omp critical
				sd_set.insert(h->helio_index);
				if (gt_sd_helio_index.count(h->helio_index)) {
#pragma omp critical
					++sd_ac;
				}
			}
			if (checkBoundingBox(cur_h->helio_pos, h->helio_pos, bk_dir, diameter)) {
#pragma omp critical
				bk_set.insert(h->helio_index);
				if (gt_bk_helio_index.count(h->helio_index)) {
#pragma omp critical
					++bk_ac;
				}
			}
		}
	}

	saveTestRes("boundingBox", sd_ac, bk_ac, sd_set, bk_set);
}


//
// [bounding box] bounding box辅助函数，判断是否有相交的可能性
//
bool SdBkCalcTest::checkBoundingBox(Vector3d& Hloc, Vector3d& HIloc, Vector3d& dir, double diameter) {
	Vector3d dist = HIloc - Hloc;
	double proj = dist.dot(dir);
	if (proj < Epsilon) return false;
	if (sqrt(pow(dist.norm(), 2) - pow(proj, 2) > diameter)) return false;
	return true;
}


//
// [neighbor Row] 对当前定日镜仅考虑当前行和相邻两行的定日镜，由东西、南北方向判断是否有可能存在阴影遮挡
//		论文表述不够清晰，所以在使用东西南北方向判断时根据bounding box
void SdBkCalcTest::neighRowSdBk() {
	Heliostat* cur_h = solar_scene->helios[helio_index];
	set<int> sdbk_set;
	int sdbk_ac = 0;
	double diameter = sqrt(pow(cur_h->helio_size.x(), 2) + pow(cur_h->helio_size.z(), 2));
	int start = 0;
	int end = 0;
	getStartEndIndex(cur_h, start, end);
#pragma omp parallel for
	for (; start <= end; ++start) {
		auto h = solar_scene->helios[start];
		if (checkEffectRegion(Vector3d(1, 0, 0), cur_h->helio_pos, h->helio_pos, diameter) ||
			checkEffectRegion(Vector3d(-1, 0, 0), cur_h->helio_pos, h->helio_pos, diameter) ||
			checkEffectRegion(Vector3d(0, 0, 1), cur_h->helio_pos, h->helio_pos, diameter) ||
			checkEffectRegion(Vector3d(0, 0, -1), cur_h->helio_pos, h->helio_pos, diameter)) {
#pragma omp critical
			sdbk_set.insert(h->helio_index);
			if (gt_sd_helio_index.count(h->helio_index) || gt_bk_helio_index.count(h->helio_index)) {
#pragma omp critical
				++sdbk_ac;
			}
		}
		
	}
	
	fstream outFile(save_path + "/neighborRow_pr.txt", ios_base::app || ios_base::out);
	outFile << sdbk_ac << ' ' << sdbk_set.size() << endl;
	outFile.close();

	outFile.open(save_path + "/neighborRow_sdbk_index.txt", ios_base::app || ios_base::out);
	for (auto&n : sdbk_set)
		outFile << n << endl;
	outFile.close();

}


//
// [neighbor row] neighbor row的辅助函数，用于确定当前定日镜当前及相邻两行的定日镜起始和终止坐标
//
void SdBkCalcTest::getStartEndIndex(Heliostat* helio, int& start, int& end) {
	vector<Heliostat*>& helios = solar_scene->helios;
	switch (solar_scene->layouts[0]->layout_type)
	{
	case RectLayoutType:
	case CrossRectLayoutType:
		double cur_row = helio->helio_pos.z();
		int cnt = 0;
		start = helio->helio_index;
		double gap = 2 * solar_scene->layouts[0]->helio_interval.z();
		while (start >= 0 && helios[start]->helio_pos.z() >= cur_row - gap) --start;
		if (start != helio->helio_index)
			++start;
		end = helio->helio_index;
		while (end < helios.size() && helios[end]->helio_pos.z() <= cur_row + gap) ++end;
		if (end != helio->helio_index)
			--end;
		break;
	case FermatLayoutType:
		double cur_row = helio->helio_poly_pos.z();
		vector<MatrixXd*>& helio_index_store = solar_scene->layouts[0]->helio_index_store;
		int cur_region = 0;
		int sum = 0;
		while (helio->helio_index > sum - 1) {
			int size = helio_index_store[cur_region]->rows() *helio_index_store[cur_region]->cols();
			sum += size;
			++cur_region;
		}
		sum -= helio_index_store[cur_region]->rows() *helio_index_store[cur_region]->cols();
		--cur_region;
		start = helio->helio_index - sum + 1;
		int row = start / helio_index_store[cur_region]->cols();
		if (row < 2 && cur_region == 0) start = 0;
		else {
			row -= 2;
			if (row < 0) {
				--cur_region;
				row %= helio_index_store[cur_region]->rows();
			}
			start = (*helio_index_store[cur_region])(row, 0);
		}
		if (row > helio_index_store[cur_region]->rows() - 3 && cur_region == helio_index_store.size() - 1) {
			end = helios.size() - 1;

		}
		else {
			row += 2;
			if (row > helio_index_store[cur_region]->rows() - 1) {
				row %= helio_index_store[cur_region]->rows();
				++cur_region;
			}
			int col = helio_index_store[cur_region]->cols() - 1;
			end = (*helio_index_store[cur_region])(row, col);
		}
		//if (row < 2) {
		//	if (cur_region == 0)
		//		start = 0;
		//	else {
		//		--cur_region;
		//		row = helio_index_store[cur_region]->rows() - 2 + row;
		//		start = (*helio_index_store[cur_region])(row, 0);
		//	}
		//}
		//else
		//	start = (*helio_index_store[cur_region])(row - 2, 0);

		//if (row > helio_index_store[cur_region]->rows() - 3) {
		//	if (cur_region == helio_index_store.size() - 1)
		//		end = helios.size() - 1;
		//	else {
		//		row = (row + 2)% helio_index_store[cur_region]->rows();
		//		++cur_region;
		//		end = (*helio_index_store[cur_region])(row, helio_index_store[cur_region]->cols() - 1);
		//	}
		//}
		//else
		//	end = (*helio_index_store)
		break;
	case RadialLayoutType:

		break;
	default:
		break;
	}
}


//
// [neighbor row] neighbor row的辅助函数，通过东西和南北的方向确定相关定日镜
//
bool SdBkCalcTest::checkEffectRegion(Vector3d dir, Vector3d& Hloc, Vector3d& HIloc, double diameter) {
	Vector3d dist = HIloc - Hloc;
	double proj = dist.dot(dir);
	if (sqrt(pow(dist.norm(), 2) - pow(proj, 2) > diameter)) return false;
	return true;
}

//
// [improved neighbor row] 在neighbor row操作之后，再进行一次剔除操作
//
void SdBkCalcTest::improvedNeighSdBk() {
	Heliostat* cur_h = solar_scene->helios[helio_index];
	set<int> sd_set, bk_set;
	int sd_ac = 0, bk_ac = 0;
	double diameter = sqrt(pow(cur_h->helio_size.x(), 2) + pow(cur_h->helio_size.z(), 2));
	int start = 0;
	int end = 0;
	getStartEndIndex(cur_h, start, end);
#pragma omp parallel for
	for (; start <= end; ++start) {
		auto h = solar_scene->helios[start];
		if (checkEffectRegion(Vector3d(1, 0, 0), cur_h->helio_pos, h->helio_pos, diameter) ||
			checkEffectRegion(Vector3d(-1, 0, 0), cur_h->helio_pos, h->helio_pos, diameter) ||
			checkEffectRegion(Vector3d(0, 0, 1), cur_h->helio_pos, h->helio_pos, diameter) ||
			checkEffectRegion(Vector3d(0, 0, -1), cur_h->helio_pos, h->helio_pos, diameter)) {
			if (checkCenterDist(cur_h, h, sd_dir)) {
				if (gt_sd_helio_index.count(h->helio_index)) {
#pragma omp critical
					++sd_ac;
				}
#pragma omp critical
				sd_set.insert(h->helio_index);
			}
			if (checkCenterDist(cur_h, h, bk_dir)) {
				if (gt_bk_helio_index.count(h->helio_index)) {
#pragma omp critical
					++bk_ac;
				}
#pragma omp critical
				bk_set.insert(h->helio_index);
			}
		}

	}

	saveTestRes("iNeighRow", sd_ac, bk_ac, sd_set, bk_set);
}


//
// [improved neighbor row] improved neighbor row辅助函数，用于计算投影坐标点，并剔除无关定日镜
//
bool SdBkCalcTest::checkCenterDist(Heliostat* H, Heliostat* HI, Vector3d& dir) {
	vector<Vector3d> helio_v, local_v, tmp_v(4);
	vector<Vector2d> project_v(4);
	double t;
	helio_v = H->vertex;
	Vector3d reverse_dir = Vector3d(-dir.x(), -dir.y(), -dir.z());

	for (int i = 0; i < helio_v.size(); i++)
		local_v.push_back(GeometryFunc::mulMatrix(helio_v[i], H->world2localM));


	vector<Vector3d> pro(4);
	helio_v = HI->vertex;
	int cnt = 0;
	for (int i = 0; i < helio_v.size(); i++) {
		t = calcIntersectionPoint(helio_v[i], reverse_dir, H->vertex[0], H->vertex[1], H->vertex[2]);
		if (t < Epsilon)
			return false;
		pro[i].x() = helio_v[i].x() + t*reverse_dir.x();
		pro[i].y() = helio_v[i].y() + t*reverse_dir.y();
		pro[i].z() = helio_v[i].z() + t*reverse_dir.z();
	}
	double phi_x = INT_MIN;
	double phi_y = INT_MIN;
	vector<Vector3d> local_proj_v;
	if (t >= Epsilon) {
		Path clip;
		for (auto v : pro) {
			local_proj_v.push_back(GeometryFunc::mulMatrix(v, H->world2localM));
			phi_x = max(local_proj_v.back().x(), phi_x);
			phi_y = max(local_proj_v.back().z(), phi_y);
		}
	}
	double h_l = H->helio_size.x();
	double h_w = H->helio_size.z();
	phi_x = 2 * phi_x / h_l;
	phi_y = 2 * phi_y / h_w;
	for (auto& v : local_proj_v) {
		if (abs(v.x()) > (phi_x + 1) / 2 * h_l || abs(v.z()) > (phi_y + 1) / 2 * h_w) return false;
	}
	return true;
		
}


//
// [3DDDA] 使用3DDDA计算光线与周围定日镜的关系
//		使用3DDDA确定相关定日镜所在网格，使用bounding sphere剔除无关定日镜
void SdBkCalcTest::use3dddaSdBk() {
	set<int> sd_set, bk_set;
	int sd_ac = 0;
	int bk_ac = 0;
	
	checkEstimateHelio(sd_dir, sd_set, sd_ac);
	checkEstimateHelio(bk_dir, bk_set, bk_ac);


	saveTestRes("3DDDA_bk", sd_ac, bk_ac, sd_set, bk_set);
}

//
// [3DDDA] 3DDDA辅助函数
//
void SdBkCalcTest::checkEstimateHelio(Vector3d& dir, set<int>& helio_set, int& ac) {
	Heliostat* helio = solar_scene->helios[helio_index];
	set<vector<int>> estimate_grid;
	calcIntersection3DDDA(helio, dir, estimate_grid);
	
	double diameter = sqrt(pow(helio->helio_size.x(), 2) + pow(helio->helio_size.z(), 2));

#pragma omp parallel for
	for (auto& iter = estimate_grid.begin(); iter != estimate_grid.end(); ++iter) {
		for (auto& rela_helio : solar_scene->layouts[0]->helio_layout[(*iter)[0]][(*iter)[1]]) {
			if (rela_helio == helio) continue;
			if (checkBoundingBox(helio->helio_pos, rela_helio->helio_pos, sd_dir, diameter)) {
				if (gt_sd_helio_index.count(helio->helio_index)) ++ac;
				helio_set.insert(rela_helio->helio_index);
			}
		}
	}
}

void SdBkCalcTest::saveTestRes(const string& file_name, const int sd_ac, const int bk_ac, const set<int>& sd_set, const set<int>& bk_set) {
	fstream outFile(save_path + "/" + file_name + "_pr.txt", ios_base::app || ios_base::out);
	outFile << sd_ac << ' ' << sd_set.size() << ' ' << bk_ac << ' ' << bk_set.size() << endl;
	outFile.close();

	outFile.open(save_path + "/" + file_name + "_sd_index.txt", ios_base::app || ios_base::out);
	for (auto&n : sd_set)
		outFile << n << endl;
	outFile.close();

	outFile.open(save_path + "/" +file_name + "_bk_index.txt", ios_base::app || ios_base::out);
	for (auto&n : bk_set)
		outFile << n << endl;
	outFile.close();
}