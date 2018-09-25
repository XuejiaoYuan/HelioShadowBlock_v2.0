
#include "./sdbkCalculator/sdbkCalculator.h"
#include "./DataStructure/SunRay.h"
#include "./LSPIA/LSPIA.h"

int main(int argc, char** argv) {

	// 计算阴影与遮挡参数
	// scene_file 镜场布置文件
	// sunray_file 太阳参数文件
	// options: -a_c 计算每个时刻下所有镜子的阴影遮挡结果(clipper);
	//			-a_r 计算每个时刻下所有镜子的阴影遮挡结果（raytracing）;
	//			-s_l 计算每个时刻下所有镜子的阴影遮挡结果（LSPIA）
	// row col 镜场行列数
	// outfile_options: -f 存储计算阴影遮挡的结果(clipper)
	//					-o 不存储结果
	// save_path 输出文件的存储地址
	if (argc != 8) {
		cout << argc << endl;
		cout << "Usage: [scene_file] [sunray_file] [options] [row] [col] [outfile_options] [save_path]" << endl;
		return -1;
	}


	string scene_filepath = string(argv[1]);
	string sunray_filepath = string(argv[2]);
	string options = string(argv[3]);
	int row = atoi(argv[4]);
	int col = atoi(argv[5]);
	string outfile_options = string(argv[6]);
	string save_path = string(argv[7]);

	SunRay sunray;
	Vector3f sunray_dir = sunray.calcSunRay(sunray_filepath);
	MatrixXf a;

	SolarScene* solar_scene = new SolarScene();
	bool flag = solar_scene->initSolarScene(scene_filepath, sunray_dir);
	if (!flag) {
		cout << "Solar scene initialize failed!" << endl;
		return -1;
	}

	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	if (options == "-s_l")
		sdbk_calc->sample_calc_preprocess(100, 100, true, true);

	vector<int> time_param(4);
	vector<MatrixXf*> gt_res;
	vector<MatrixXf*> sample_sd_bk_res;
	auto start = std::chrono::high_resolution_clock::now();
	//for (int month = 1; month <= 12; month++) {
	//	for (int day = 1; day < 29; day += 4) {
	//		for (int hour = 8; hour < 13; hour++) {
	//			for (int min = 0; min < 15; min += 15) {
	//				time_param[0] = month;
	//				time_param[1] = day;
	//				time_param[2] = hour;
	//				time_param[3] = min;
	//				sunray_dir = sunray.changeSunRay(time_param);
	//				solar_scene->changeSolarScene(sunray_dir);
	//				
	//				if (options == "-a_c") {
	//					gt_res.push_back(sdbk_calc->calcShadowBlock());
	//					if(outfile_options=="-f")
	//						sdbk_calc->save_clipper_res(save_path, month, day, hour, min);
	//				}
	//				else if (options == "-a_r")
	//					gt_res.push_back(sdbk_calc->calcShadowBlock());
	//				else {
	//					gt_res.push_back(sdbk_calc->calcShadowBlock());
	//					sample_sd_bk_res.push_back(sdbk_calc->calcSampleShadowBlock());
	//				}
	//			}
	//		}
	//	}
	//}



	time_param[0] = 1;
	time_param[1] = 1;
	time_param[2] = 8;
	time_param[3] = 0;
	sunray_dir = sunray.changeSunRay(time_param);
	solar_scene->changeSolarScene(sunray_dir);

	if(options=="-s_l") {
		gt_res.push_back(sdbk_calc->calcShadowBlock());
		sample_sd_bk_res.push_back(sdbk_calc->calcSampleShadowBlock());
	}

	for (int i = 0; i < sdbk_calc->sample_field_data[0]->rows(); i++)
		cout << (*sdbk_calc->sample_field_data[1])(i, 0) << endl;

	vector<vector<MatrixXf*>> fitting_sd_bk_res;
	if (options == "-s_l") {
		LSPIA lspia;
		lspia.set_datas(sdbk_calc->field_data, sdbk_calc->sample_field_data, sample_sd_bk_res);

		vector<int> ctrl_nums = { 76,76 };
		fitting_sd_bk_res = lspia.LSPIA_surface(ctrl_nums, 0.9);
	}
	//int sample_row_num = 30;
	//int sample_col_num = 30;

	//for (int i = 0; i < sample_row_num / 2 - 1; i++) {
	//	int f_i = int(row*i / sample_row_num);
	//	for (int k = 0; k < 2; k++) {
	//		for (int j = 0; j < sample_col_num - 1; j++) {
	//			int f_j = int(col*j / sample_col_num);
	//			fstream outFile(save_path + "sample/h_fi" + to_string(f_i) + "_fj_" + to_string(f_j) + ".txt", ios_base::out);
	//			for (int t = 0; t < gt_res.size(); t++) {
	//				outFile << (*gt_res[t])(2 * f_i + k, f_j) << endl;
	//			}
	//		}
	//	}
	//}


	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	std::cout << "Total clipper time: " << time << "s." << endl;

	return 1;
}