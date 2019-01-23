
#include "./sdbkCalculator/sdbkCalculator.h"
#include "./DataStructure/SunRay.h"
#include "./LSPIA/LSPIA.h"
#include "./GaussLegendre/GaussLegendre.h"
#include "./DataStructure/FieldSegment.h"
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
	Vector3d sunray_dir = sunray.calcSunRay(sunray_filepath);
	MatrixXd a;

	SolarScene* solar_scene = new SolarScene();
	solar_scene->initFieldParam(scene_filepath);

	// Gauss-Legendre参数预计算
	int N = solar_scene->recvs[0]->recv_size.x() / RECEIVER_SLICE;
	int M = solar_scene->recvs[0]->recv_size.z() / RECEIVER_SLICE;
	GaussLegendre* gl = new GaussLegendre();
	N = 8;
	M = 8;
	gl->calcNodeWeight(N, M);

	// TODO 镜场参数优化
	int cnt;
	vector<vector<double>*> field_args;
	switch (solar_scene->layouts[0]->layout_type)
	{
	case CrossRectFieldType:
	case RectLayoutType:
		field_args.push_back(new vector<double>{10, 10, 10});			// 定日镜间隔
		field_args.push_back(new vector<double>{ -80 });					// 第一行定日镜与接收器之间的距离
		field_args.push_back(new vector<double>{ 300 });					// 定日镜行数
		field_args.push_back(new vector<double>{ 300 });					// 定日镜列数
		break;
	case FermatLayoutType:
		field_args.push_back(new vector<double>{ double(0.1) });			// 定日镜包围盒安全距离
		field_args.push_back(new vector<double>{ double( 114.9417 )});	// 第一个同心圆与接收器之间的距离(Campo: 114.9417)
		field_args.push_back(new vector<double>{ 0.866f });				// 第一个同心环中定日镜的分布间隔
		field_args.push_back(new vector<double>{ 1.4f });				// 第二个同心环中定日镜的分布间隔
		field_args.push_back(new vector<double>{ 2.0f});					// 第三个同心环中定日镜的分布间隔
		break;
	default:
		break;
	}
	solar_scene->adjustFieldParam(field_args);
	//FieldSegmentCreator field_seg_creator;
	//FieldSegment* field_seg = field_seg_creator.getFieldSegment(solar_scene->layouts[0]->layout_type);
	//field_seg->initFieldMatrix(solar_scene);
	//field_seg->updateSample(10, 10);
	FieldSegment* field_seg = new FieldSegment(solar_scene);
	field_seg->setSegmentParam(10, 10, 100, 100);
	field_seg->initFieldSegment();

	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene, gl);


	// 最小二乘拟合数据
	//vector<vector<int>> daily = {
	//	{ 3,21 },		// 春分
	//	{ 6,22 },		// 夏至
	//	{ 9,23 },		// 秋分
	//	{ 12,22 }		// 冬至
	//};
	//vector<int> sample_hour = { 8, 10, 12, 14, 16 };

	//for (int helio_index = 0; helio_index < solar_scene->helios.size(); helio_index+=100) {
	//	fstream outFile(save_path + "/single_helio_" + to_string(helio_index) + ".txt", ios_base::out);
	//	//fstream outFile(save_path + "/sunray_altitude.txt", ios_base::out);
	//	auto start = std::chrono::high_resolution_clock::now();
	//	for (int i = 0; i < daily.size(); i++) {
	//		for (int j = 0; j < sample_hour.size(); j++) {
	//			vector<int> t = { daily[i][0], daily[i][1], sample_hour[j], 0 };
	//			sunray_dir = sunray.changeSunRay(t);
	//			solar_scene->changeSolarScene(sunray_dir);
	//			double visi = 1- sdbk_calc->calcSingleShadowBlock(helio_index);
	//			double flux_sum = sdbk_calc->calcSingleFluxSum(helio_index, sunray.current_DNI);
	//			outFile << setprecision(12) << visi << ' ' << flux_sum << ' ' << visi*flux_sum << endl;
	//		}
	//		//outFile << endl;
	//	}
	//	outFile.close();
	//	//auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	//	//auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	//	//std::cout << "Total clipper time: " << time << "s." << endl;

	//	outFile.open(save_path + "/total_single_helio_" + to_string(helio_index) + ".txt", ios_base::out);
	//	vector<int> time_param(4);
	//	double sum = 0;
	//	double cnt = 0;
	//	double dni_sum = 0;
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
	//					double visi = 1 - sdbk_calc->calcSingleShadowBlock(helio_index);
	//					double flux_sum = sdbk_calc->calcSingleFluxSum(helio_index, sunray.current_DNI);
	//					outFile << setprecision(12) << visi << ' ' << flux_sum << ' ' << visi*flux_sum << endl;
	//					cnt++;
	//				}
	//			}
	//			//outFile << endl;
	//		}
	//	}
	//	outFile.close();
	//}

	vector<int> time_param(4, 0);
	//time_param[0] = 6;
	//time_param[1] = 22;
	//time_param[2] = 12;
	//time_param[3] = 0;
	//sunray_dir = sunray.changeSunRay(time_param);
	//solar_scene->changeSolarScene(sunray_dir);

	//sdbk_calc->calcSingleFluxSum(0, sunray.current_DNI);
	//return 0;

	vector<MatrixXd*> gt_res;
	double total_t = 0;
	cnt = 0;
	vector<double> calc_t;
	//fstream sunray_file("sunray_dir.txt", ios_base::out | ios_base::app);
	//for (int month = 1; month <= 6; month++) {
	//	for (int day = 1; day < 29; day += 4) {
	//		for (int hour = 8; hour < 13; hour++) {
	//			for (int min = 0; min < 15; min += 15) {
	//				time_param[0] = month;
	//				time_param[1] = day;
	//				time_param[2] = hour;
	//				time_param[3] = min;

	//				sunray_dir = sunray.changeSunRay(time_param);
	//				sunray_file << sunray_dir.x() << ' ' << sunray_dir.z() << endl;
	//			}
	//		}
	//	}
	//}
	//sunray_file.close();

	//for (int month = 1; month <= 6; month++) {
	//	for (int day = 1; day < 29; day += 4) {
	//		for (int hour = 8; hour < 13; hour++) {
	//			for (int min = 0; min < 15; min += 15) {
	//				cnt++;
	//				cout << "cnt: " << cnt << endl;
	//				auto start = std::chrono::high_resolution_clock::now();

	//				time_param[0] = month;
	//				time_param[1] = day;
	//				time_param[2] = hour;
	//				time_param[3] = min;
	//				
	//				sunray_dir = sunray.changeSunRay(time_param);
	//				solar_scene->changeSolarScene(sunray_dir);
	//				
	//				if (options == "-a_c") {
	//					//gt_res.push_back(sdbk_calc->calcShadowBlock());
	//					sdbk_calc->calcShadowBlock(sunray.current_DNI);
	//					if (outfile_options == "-f")
	//						sdbk_calc->save_clipper_res(save_path, month, day, hour, min);
	//				}
	//				else if (options == "-a_r")
	//					//gt_res.push_back(sdbk_calc->calcShadowBlock());
	//					sdbk_calc->calcShadowBlock(sunray.current_DNI);
	//				else
	//					gt_res.push_back(sdbk_calc->calcSampleShadowBlock());
	//				auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	//				auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	//				total_t += time;
	//				calc_t.push_back(time);
	//				std::cout << "Change and Calculate time: " << time << "s." << endl;

	//			}
	//		}
	//	}
	//}
	//cout << "average time: " << total_t / cnt << " time cnt: " << cnt << endl;
	//fstream outFile("calc_t.txt", ios_base::out);
	//for (auto& t : calc_t) {
	//	outFile << t << endl;
	//}
	//outFile.close();
	////
	//// test
	////

	time_param[0] = 1;
	time_param[1] = 1;
	time_param[2] = 12;
	time_param[3] = 0;

	auto start = std::chrono::high_resolution_clock::now();
	sunray_dir = sunray.changeSunRay(time_param);
	solar_scene->changeSolarScene(sunray_dir);
	sdbk_calc->calcShadowBlock(sunray.current_DNI);

	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	std::cout << "Total clipper time: " << time << "s." << endl;

	start = std::chrono::high_resolution_clock::now();

	//if (options == "-s_l")
	//	sdbk_calc->sample_calc_preprocess(100, 100, true, true);

	//vector<vector<MatrixXd*>> sample_sd_bk_res;
	//sunray_dir = sunray.changeSunRay(time_param);
	//solar_scene->changeSolarScene(sunray_dir);
	if (options == "-s_l") {
		field_seg->even_sample_res.push_back(sdbk_calc->calcSampleShadowBlock(field_seg->even_sample_field_index, sunray.current_DNI));
		field_seg->odd_sample_res.push_back(sdbk_calc->calcSampleShadowBlock(field_seg->odd_sample_field_index, sunray.current_DNI));
		sdbk_calc->calcExcludeShadowBlock(sunray.current_DNI);
	}

	vector<vector<MatrixXd*>> fitting_sd_bk_res;
	if (options == "-s_l") {
		LSPIA lspia;
		//lspia.set_datas(sdbk_calc->field_data, sdbk_calc->sample_field_data, sample_sd_bk_res);

		vector<int> ctrl_nums = { 50, 50 };
		//fitting_sd_bk_res = lspia.LSPIA_surface(ctrl_nums, 1.7);		// sd bk: 76,76,0.9

		//lspia.checkFittingData(sdbk_calc->solar_scene->helios, sdbk_calc->field_index, fitting_sd_bk_res);
		
		lspia.setPreDatas(field_seg, ctrl_nums, 0.1);
		lspia.LSPIA_surface();
	}

	elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	std::cout << "LSPIA calculation time: " << time << "s." << endl;

	return 1;
}