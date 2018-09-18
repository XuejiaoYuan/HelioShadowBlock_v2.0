
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
	vector<vector<MatrixXf*>> gt_res;
	vector<vector<MatrixXf*>> sample_res;
	auto start = std::chrono::high_resolution_clock::now();
	for (int month = 1; month <= 12; month++) {
		for (int day = 1; day < 29; day += 4) {
			for (int hour = 8; hour < 18; hour++) {
				for (int min = 0; min < 60; min += 15) {
					time_param[0] = month;
					time_param[1] = day;
					time_param[2] = hour;
					time_param[3] = min;
					sunray_dir = sunray.changeSunRay(time_param);
					solar_scene->changeSolarScene(sunray_dir);
					if (options == "-a_c") {
						gt_res.push_back(sdbk_calc->calcShadowBlock());
						fstream outFile(save_path + "/gt_clipper_sd&bk_m" + to_string(month) + "_d" + to_string(day) + "_h" 
							+ to_string(hour) + "_min" + to_string(min) + ".txt", ios_base::out);
						outFile << row << ' ' << col << endl;
						for(int i)
					}
					else if (options == "-a_r")
						gt_res.push_back(sdbk_calc->calcShadowBlock());
					else
						sample_res.push_back(sdbk_calc->calcSampleShadowBlock());
				}
			}
		}
	}

	//
	// test
	//
	time_param[0] = 1;
	time_param[1] = 1;
	time_param[2] = 8;
	time_param[3] = 0;
	sunray_dir = sunray.changeSunRay(time_param);
	solar_scene->changeSolarScene(sunray_dir);
	if (options == "-s_l")
		sample_res.push_back(sdbk_calc->calcSampleShadowBlock());

	if (options == "-s_l") {
		vector<int> ctrl_num = { 76, 76 };
		float miu = 1e-4;
		LSPIA* lspia = new LSPIA();
		lspia->set_datas(sdbk_calc->field_data, sdbk_calc->sample_field_data, sample_res);
		lspia->LSPIA_surface(ctrl_num, miu);
	}

	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	std::cout << "Total clipper time: " << time << "s." << endl;

	return 1;
}