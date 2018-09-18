float SolarScene::fhelioClipper(Heliostat * helio, const float3 & dir, const set<vector<int>>& estimate_grid)
{
	vector<float3> helio_v, local_v;
	float2 project_v;
	float t;
	//Paths subj(1), clips;
	PolyPList lightPart;
	Point p[4];
	helio_v = helio->vertex;
	float3 reverse_dir = make_float3(-dir.x(), -dir.y(), -dir.z());

	for (int i = 0; i < helio_v.size(); i++) {
		local_v.push_back(global_func::mulMatrix(helio_v[i], helio->world2localM));
		p[i] = Point(local_v[i].x(), local_v[i].z());
	}
	Poly origPoly = Poly(p[0]);
	for (int i = 1; i < 4; i++)
		origPoly.add(p[i]);
	lightPart.add(&origPoly);

	float total_area = origPoly.area();

	if (total_area == 0) {
		cout << "Project surface is 0!" << endl;
		return 0;
	}

	float3 local_dir = global_func::mulMatrix(reverse_dir, helio->world2localM, false);
	local_dir = normalize(local_dir);
	for (auto iter = estimate_grid.begin(); iter != estimate_grid.end(); iter++) {
		for (auto&relative_helio : layouts[0]->helio_layout[(*iter)[0]][(*iter)[1]]) {
			if (relative_helio == helio)
				continue;
			helio_v = relative_helio->vertex;
			Point project_v[4];
			int index = 0;
			for (index = 0; index < helio_v.size(); index++) {
				local_v[index] = global_func::mulMatrix(helio_v[index], helio->world2localM);
				t = -local_v[index].y() / local_dir.y();
				if (t < Epsilon)
					break;
				project_v[index].x() = local_v[index].x() + t*local_dir.x();
				project_v[index].y() = local_v[index].z() + t*local_dir.z();
			}
			if (index == helio_v.size()) {
				Poly occludeShape = Poly(project_v[0]);
				for (int i = 1; i < 4; i++)
					occludeShape.add(project_v[i]);

				PolyPList e_min_d, d_min_e, d_and_e;
				PolyPListIter iter(lightPart);
				while (iter())
					clip_poly(*iter.val(), occludeShape, e_min_d, d_min_e, d_and_e);    //kernel
																						// 		clear(d_min_e);														   // 		clear(d_and_e); 
																						//清空原光亮区形状容器
				PolyPListIter iter1(lightPart);
				while (iter1())
					lightPart.del(iter1.val());

				PolyPListIter iter2(lightPart);
				while (iter2())
					lightPart.del(iter2.val());

				//光亮区形状容器添加裁剪后多边形
				PolyPListIter iter3(e_min_d);
				while (iter3())
					lightPart.add(iter3.val());
			}
		}
	}

	float sum = 0;
	PolyPListIter iter(lightPart);
	int poly_cnt = 0;
	vector<vector<float2>> clipPath;
	while (iter()) {
		Poly tmp = Poly(*iter.val());
		ConstPolyIter cIter(tmp);
		int index = 0;
		vector<float2> v;
		while (cIter()) {
			v.push_back(make_float2(cIter.point().x(), cIter.point().y()));
		}
		clipPath.push_back(v);
		sum += tmp.area();
	}

	float res = 1 - sum / total_area;

#ifdef OUTPUTRES


#endif // OUTPUTRES

	return res;

}
