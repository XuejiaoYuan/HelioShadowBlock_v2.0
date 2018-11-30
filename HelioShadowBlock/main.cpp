
#include "./sdbkCalculator/sdbkCalculator.h"
#include "./DataStructure/SunRay.h"
#include "./LSPIA/LSPIA.h"
#include "./GaussLegendre/GaussLegendre.h"
#include <iomanip>

int get_index(int r, int c) {
	int even_rows = r / 2 + r % 2;
	int odd_rows = r / 2;
	return even_rows * 300 + odd_rows * 299 + c;
}

int main(int argc, char** argv) {

	// 计算阴影与遮挡参数
	// scene_file 镜场布置文件
	// sunray_file 太阳参数文件
	// options: -a_c 计算每个时刻下所有镜子的阴影遮挡结果(clipper);
	//			-a_r 计算每个时刻下所有镜子的阴影遮挡结果（raytracing）;
	//			-s_l 计算每个时刻下所有镜子的阴影遮挡结果（LSPIA）
	// outfile_options: -f 存储计算阴影遮挡的结果(clipper)
	//					-o 不存储结果
	// save_path 输出文件的存储地址
	if (argc != 6) {
		cout << argc << endl;
		cout << "Usage: [scene_file] [sunray_file] [options] [outfile_options] [save_path]" << endl;
		return -1;
	}


	string scene_filepath = string(argv[1]);
	string sunray_filepath = string(argv[2]);
	string options = string(argv[3]);
	string outfile_options = string(argv[4]);
	string save_path = string(argv[5]);

	SunRay sunray;
	Vector3f sunray_dir = sunray.calcSunRay(sunray_filepath);
	MatrixXf a;

	SolarScene* solar_scene = new SolarScene();
	solar_scene->initFieldParam(scene_filepath);

	// Gauss-Legendre参数预计算
	int N = solar_scene->recvs[0]->recv_size.x() / RECEIVER_SLICE;
	int M = solar_scene->recvs[0]->recv_size.z() / RECEIVER_SLICE;
	GaussLegendre* gl = new GaussLegendre();
	gl->calcNodeWeight(N, M);

	// TODO 镜场参数优化
	int cnt;
	vector<vector<float>*> field_args;
	switch (solar_scene->layouts[0]->layout_type)
	{
	case CrossRectFieldType:
	case RectLayoutType:
		field_args.push_back(new vector<float>{10, 10, 10});			// 定日镜间隔
		field_args.push_back(new vector<float>{ -80 });					// 第一行定日镜与接收器之间的距离
		field_args.push_back(new vector<float>{ 300 });					// 定日镜行数
		field_args.push_back(new vector<float>{ 300 });					// 定日镜列数
		break;
	case FermatLayoutType:
		field_args.push_back(new vector<float>{ float(0.1) });			// 定日镜包围盒安全距离
		field_args.push_back(new vector<float>{ float( 114.9417 )});	// 第一个同心圆与接收器之间的距离(Campo: 114.9417)
		field_args.push_back(new vector<float>{ 0.866f });				// 第一个同心环中定日镜的分布间隔
		field_args.push_back(new vector<float>{ 1.4f });				// 第二个同心环中定日镜的分布间隔
		field_args.push_back(new vector<float>{ 2.0f});					// 第三个同心环中定日镜的分布间隔
		break;
	default:
		break;
	}
	solar_scene->adjustFieldParam(field_args);
	//bool flag = solar_scene->initSolarScene(scene_filepath, sunray_dir);
	//if (!flag) {
	//	cout << "Solar scene initialize failed!" << endl;
	//	return -1;
	//}

	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene, gl);
	//if (options == "-s_l")
	//	sdbk_calc->sample_calc_preprocess(100, 100, true, true);


	// 最小二乘拟合数据
	//vector<int> r = { 0, 74, 149, 224, 299 };
	//vector <int> odd = { 0, 74, 148, 222, 296 };
	//vector<vector<int>> daily = {
	//	{ 3,21 },		// 春分
	//	{ 6,22 },		// 夏至
	//	{ 9,23 },		// 秋分
	//	{ 12,22 }		// 冬至
	//};
	//vector<int> sample_hour = { 8, 10, 12, 14, 16 };

	//for (int helio_index = 12700; helio_index < solar_scene->helios.size(); helio_index+=50) {
	//	fstream outFile(save_path + "/single_helio_" + to_string(helio_index) + ".txt", ios_base::out);
	//	//fstream outFile(save_path + "/sunray_altitude.txt", ios_base::out);
	//	auto start = std::chrono::high_resolution_clock::now();
	//	for (int i = 0; i < daily.size(); i++) {
	//		for (int j = 0; j < sample_hour.size(); j++) {
	//			vector<int> t = { daily[i][0], daily[i][1], sample_hour[j], 0 };
	//			sunray_dir = sunray.changeSunRay(t);
	//			solar_scene->changeSolarScene(sunray_dir);
	//			sdbk_calc->calcSingleShadowBlock(helio_index);
	//			solar_scene->helios[helio_index]->calcSunHelioAngle(sunray_dir);
	//			outFile << setprecision(12) << 1 - sdbk_calc->calcSingleShadowBlock(helio_index) << ' ' 
	//				<< sunray.current_DNI <<' ' << solar_scene->helios[helio_index]->calcSunHelioAngle(sunray_dir) << ' ';
	//		}
	//		outFile << endl;
	//	}
	//	outFile.close();
	//	//auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	//	//auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	//	//std::cout << "Total clipper time: " << time << "s." << endl;

	//	outFile.open(save_path + "/total_single_helio_" + to_string(helio_index) + ".txt", ios_base::out);
	//	vector<int> time_param(4);
	//	float sum = 0;
	//	float cnt = 0;
	//	float dni_sum = 0;
	//	for (int month = 1; month <= 12; month++) {
	//		for (int day = 1; day < 29; day += 4) {
	//			for (int hour = 8; hour < 17; hour++) {
	//				for (int min = 0; min < 30; min += 30) {
	//					time_param[0] = month;
	//					time_param[1] = day;
	//					time_param[2] = hour;
	//					time_param[3] = min;
	//					sunray_dir = sunray.changeSunRay(time_param);
	//					solar_scene->changeSolarScene(sunray_dir);
	//					float visi = 1 - sdbk_calc->calcSingleShadowBlock(helio_index);
	//					float dni = sunray.current_DNI;
	//					float cos_w = solar_scene->helios[helio_index]->calcSunHelioAngle(sunray_dir);
	//					sum += visi * dni * cos_w;
	//					dni_sum += sunray.current_DNI;
	//					outFile << setprecision(12) << visi << ' ' << dni << ' ' << cos_w << ' ';
	//					cnt++;
	//				}
	//			}
	//			outFile << endl;
	//		}
	//	}
	//	outFile << sum / dni_sum << endl;
	//	outFile.close();
	//}


	vector<MatrixXf*> gt_res;
	vector<MatrixXf*> sample_sd_bk_res;
	auto start = std::chrono::high_resolution_clock::now();
	vector<int> time_param(4, 0);
	for (int month = 1; month <= 12; month++) {
		for (int day = 1; day < 29; day += 4) {
			for (int hour = 8; hour < 13; hour++) {
				for (int min = 0; min < 15; min += 15) {
					time_param[0] = month;
					time_param[1] = day;
					time_param[2] = hour;
					time_param[3] = min;
					sunray_dir = sunray.changeSunRay(time_param);
					solar_scene->changeSolarScene(sunray_dir);

					if (options == "-a_c") {
						//gt_res.push_back(sdbk_calc->calcShadowBlock());
						sdbk_calc->calcShadowBlock(sunray.current_DNI);
						if (outfile_options == "-f")
							sdbk_calc->save_clipper_res(save_path, month, day, hour, min);
					}
					else if (options == "-a_r")
						//gt_res.push_back(sdbk_calc->calcShadowBlock());
						sdbk_calc->calcShadowBlock(sunray.current_DNI);
					else
						gt_res.push_back(sdbk_calc->calcSampleShadowBlock());
				}
			}
		}
	}
	////
	//// test
	////
	//time_param[0] = 1;
	//time_param[1] = 1;
	//time_param[2] = 8;
	//time_param[3] = 0;
	//sunray_dir = sunray.changeSunRay(time_param);
	//solar_scene->changeSolarScene(sunray_dir);
	//if (options == "-s_l")
	//	sample_sd_bk_res.push_back(sdbk_calc->calcSampleShadowBlock());

	//for (int i = 0; i < sdbk_calc->sample_field_data[0]->rows(); i++)
	//	cout << (*sdbk_calc->sample_field_data[1])(i, 0) << endl;

	//vector<vector<MatrixXf*>> fitting_sd_bk_res;
	//if (options == "-s_l") {
	//	LSPIA lspia;
	//	lspia.set_datas(sdbk_calc->field_data, sdbk_calc->sample_field_data, sample_sd_bk_res);

	//	vector<int> ctrl_nums = { 76,76 };
	//	fitting_sd_bk_res = lspia.LSPIA_surface(ctrl_nums, 0.9);
	//}

	//auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	//auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	//std::cout << "Total clipper time: " << time << "s." << endl;

	return 1;
}