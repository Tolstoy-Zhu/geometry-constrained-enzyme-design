"""Additional controls inside the original StructureMine workbench."""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
import streamlit as st
from .families import profiles,get_profile,reference
from .store import Store,state_root,index_history,uid
from .execution import submit,server_profiles,ServerProfile,save_server,RemoteExecutor
from .services import import_records,export_assays,reconcile_task,retry_task,cancel_task,register_model,register_task_model

def context():
    return {"project_id":st.session_state.get("aqua_project","default"),
            "family_id":st.session_state.get("aqua_family","bur01"),
            "executor":st.session_state.get("aqua_executor","local"),
            "reaction_id":st.session_state.get("aqua_reaction") or None}

def project_sidebar():
    store=Store();rows=store.projects();names={p["id"]:p["name"] for p in rows}
    if "_aqua_project_seen" not in st.session_state and st.query_params.get("run"):
        try:
            linked=store.task(st.query_params["run"])["project_id"]
            if linked in names:st.session_state["aqua_project"]=linked
        except ValueError:pass  # Historical filesystem report URLs remain valid in the default project.
    with st.sidebar:
        st.subheader("项目与计算位置")
        project=st.selectbox("当前项目",list(names),format_func=names.__getitem__,key="aqua_project")
        selected=next(p for p in rows if p["id"]==project)
        if st.session_state.get("_aqua_project_seen")!=project:
            previous=st.session_state.get("_aqua_project_seen")
            drafts=st.session_state.setdefault("_aqua_project_drafts",{})
            fields=("web_fasta","web_chem","web_parent","web_source","web_batches","web_report_choice",
                    "kinetic_fasta","kinetic_smiles","kinetic_enzyme_records","kinetic_model","kinetic_report")
            if previous:
                drafts[previous]={k:st.session_state[k] for k in fields if k in st.session_state}
                if "run" in st.query_params:del st.query_params["run"]
                for key in fields:st.session_state.pop(key,None)
                st.session_state.update(drafts.get(project,{}))
                st.session_state.pop("_web_select_job",None)
            st.session_state["aqua_family"]=selected["family_id"]
            st.session_state["aqua_executor"]=selected["executor"]
            st.session_state["_aqua_project_seen"]=project
        families={p.family_id:p.name for p in profiles()}
        family=st.selectbox("酶家族",list(families),format_func=families.__getitem__,key="aqua_family")
        servers=server_profiles();executors=["local",*servers]
        if st.session_state.get("aqua_executor") not in executors:st.session_state["aqua_executor"]="local"
        st.selectbox("计算位置",executors,format_func=lambda key:"本机" if key=="local" else key,key="aqua_executor")
        reactions={r["reaction_id"]:r["name"] for r in store.records("reaction",project)}
        if st.session_state.get("aqua_reaction") not in reactions:st.session_state["aqua_reaction"]=""
        st.selectbox("关联反应条件",["",*reactions],format_func=lambda k:reactions.get(k,"未指定"),key="aqua_reaction")
        p=get_profile(family)
        role_names={"nucleophile":"亲核残基","acid":"酸性残基","base":"催化碱","acid_base":"酸碱催化残基"}
        st.caption("催化位点："+ " · ".join(f"{s.amino_acid}{s.residue_number}（{role_names.get(s.role,s.role)}）" for s in p.sites))
        if not p.wase_screen:st.caption("此家族单独检查结构与催化位点，不计算 Bur01 的 WASE 分数。")
        with st.expander("新建项目"):
            name=st.text_input("项目名称",key="aqua_new_project")
            if st.button("创建项目",key="aqua_create_project"):
                identifier=store.create_project(name,family,next(iter(servers),"local"))
                st.session_state["_aqua_select_project"]=identifier
                st.success("项目已创建，可在上方选择。")
    return p

def _do(operation):
    try:return operation()
    except (ValueError,OSError,RuntimeError,KeyError,TypeError) as exc:
        st.error(str(exc));return None

def _submit(kind,inputs):
    identifier=submit(kind,inputs,**context())
    st.success("任务已保存："+identifier+"。进度与下载见“项目与任务”。")
    return identifier

@st.fragment(run_every="10s")
def project_page():
    st.subheader("项目与任务")
    st.caption("任务和输入保存到本机索引；远端计算不依赖浏览器保持打开。")
    ctx=context();store=Store()
    from .polling import schedule_sync
    schedule_sync(ctx["project_id"],store)
    enzymes=store.records("enzyme",ctx["project_id"])
    if enzymes:
        with st.expander("酶与表达构建",expanded=False):
            st.dataframe(pd.DataFrame([{"名称":r["name"],"编号":r["enzyme_id"],"构建":r["construct_id"],
                "长度":len(r["sequence"]),"亲本":r.get("parent_id"),"序列摘要":r["sequence_sha256"]} for r in enzymes]),
                hide_index=True,width="stretch")
            selected=st.selectbox("查看构建记录",[r["enzyme_id"] for r in enzymes])
            record=next(r for r in enzymes if r["enzyme_id"]==selected)
            st.json(record)
            st.download_button("下载此构建序列",">"+record["enzyme_id"]+"\n"+record["sequence"]+"\n",
                               file_name=record["enzyme_id"]+".fasta")
    with st.expander("计算服务器配置"):
        st.caption("连接设置只保存在本机。使用已配置的 SSH 密钥；此处不保存密码。")
        with st.form("aqua_server"):
            name=st.text_input("服务器名称",value="autodl")
            host=st.text_input("SSH 地址")
            user=st.text_input("用户名",value="root")
            port=st.number_input("SSH 端口",1,65535,22)
            identity=st.text_input("本地密钥文件路径（已用 SSH 默认密钥时可留空）")
            root=st.text_input("远端工作目录",value="/root/autodl-tmp/aquaesternet")
            interpreter=st.text_input("远端环境 Python",value="/root/miniconda3/envs/aqua/bin/python")
            if st.form_submit_button("保存服务器"):
                saved=_do(lambda:save_server(ServerProfile(name=name,host=host,user=user,port=int(port),
                            identity_file=identity or None,remote_root=root,python=interpreter)))
                if host:st.caption("请使用“检查连接”确认配置和远端程序。")
        servers=server_profiles()
        if servers:
            selected=st.selectbox("检查的服务器",list(servers),key="server_check_name")
            if st.button("检查连接"):
                result=_do(lambda:RemoteExecutor(servers[selected]).check())
                if result:st.json(result)
    rows=store.tasks(ctx["project_id"])
    if rows:
        frame=pd.DataFrame([{"编号":r["id"],"任务":json.loads(r["spec"])["kind"],"状态":r["state"],
                             "计算位置":json.loads(r["spec"])["executor"],"更新时间":r["updated"],"说明":r["message"]} for r in rows])
        st.dataframe(frame,hide_index=True,width="stretch")
        task=st.selectbox("查看任务",[r["id"] for r in rows],key="aqua_task")
        row=store.task(task)
        a,b,c=st.columns(3)
        if a.button("同步状态 / 重新连接"):
            if _do(lambda:reconcile_task(task,store)) is not None:st.rerun()
        if b.button("取消此任务",disabled=row["state"] in {"completed","cancelled"}):
            _do(lambda:cancel_task(task,ctx["project_id"],store));st.rerun()
        if c.button("重试已停止任务",disabled=row["state"] not in {"failed","interrupted","cancelled"}):
            _do(lambda:retry_task(task,store));st.rerun()
        directory=store.path.parent/"tasks"/task
        log=directory/"process.log"
        if log.exists():st.code(log.read_text(encoding="utf-8",errors="replace")[-18000:],language=None)
        if row.get("result"):
            result=json.loads(row["result"]);st.json(result)
            output=Path(result.get("local_output",result.get("output",directory/"output")))
            files=[p for p in output.rglob("*") if p.is_file() and p.suffix in {".json",".csv",".tsv",".html",".fasta",".xlsx",".pdb",".zip"}]
            if files:
                file=st.selectbox("下载结果文件",files,format_func=lambda p:str(p.relative_to(output)),key="aqua_download")
                if file.stat().st_size<=50*1024*1024:
                    st.download_button("下载",file.read_bytes(),file_name=file.name,key="aqua_download_button")
    else:st.info("此项目尚无后台任务。")
    with st.expander("历史结果索引"):
        from ..config import project_root
        root=st.text_input("既有仓库目录",value=str(project_root()),key="history_root")
        if st.button("建立索引"):
            count=_do(lambda:index_history(root,store,ctx["project_id"]))
            if count is not None:st.success(f"已索引 {count} 份记录，原始文件保留。")
        records=store.records("history",ctx["project_id"])
        if records:st.dataframe(pd.DataFrame(records),hide_index=True)
    with st.expander("实验记录入库"):
        st.caption("导入酶、反应和实测记录，保留方法、单位、重复编号、检测限及原始来源。")
        upload=st.file_uploader("实验记录 JSON",type=["json"],key="assay_records")
        if upload and st.button("校验并入库"):
            path=state_root()/"imports"/(uid("assays")+".json");path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(upload.getvalue())
            result=_do(lambda:import_records(path,ctx["project_id"],store))
            if result:st.success(str(result))
        if st.button("导出动力学数据表"):
            output=state_root()/"exports"/(uid("assays")+".csv")
            result=_do(lambda:export_assays(output,ctx["project_id"],store))
            if result:st.download_button("下载实测数据表",output.read_bytes(),file_name="assays.csv")
        assays=store.records("assay",ctx["project_id"])
        if assays:st.dataframe(pd.DataFrame(assays),hide_index=True,width="stretch")

def model_page():
    from .polling import schedule_sync
    from .services import kinetic_reports
    schedule_sync(context()['project_id'])
    models=Store().records('model',context()['project_id'])
    st.subheader("模型评估与训练")
    st.write("kcat 和 Km 分别训练；酯合成初速、水解初速、合成/水解比目前尚未训练。")
    st.caption("冻结分组：70% 训练、10% 验证、10% 区间校准、10% 测试。同源序列及同一亲本的记录保持同组。")
    if models:
        st.dataframe(pd.DataFrame([{'指标':r['card']['endpoint'],'模型':r['card']['model'],
            '测试样本':r['card']['metrics']['new_enzyme']['n'],
            '对数 MAE':r['card']['metrics']['new_enzyme']['MAE'],
            'Ridge 对照 MAE':r['card']['metrics']['Ridge']['MAE'],
            '发布状态':'通过对照评估' if r['card']['status']=='validated_upgrade' else '实验模型'} for r in models]),
            hide_index=True,width='stretch')
    with st.expander("构建公开数据或实验数据",expanded=not models):
        source=st.selectbox("数据来源",["catapro","csv"],format_func=lambda k:"CataPro 冻结公开数据" if k=="catapro" else "已映射的实测 CSV")
        path=st.text_input("数据文件路径",key="dataset_input")
        endpoint=st.selectbox("训练指标",["kcat","Km"])
        cluster=st.text_input("已有同源分组 TSV（可留空）")
        auto_cluster=st.checkbox("自动完成 40% 一致性、80% 覆盖度分组",value=True)
        if st.button("构建数据集"):_do(lambda:_submit("dataset",{"input":path,"source":source,
            "endpoint":endpoint,"cluster_file":cluster or None,"mmseqs":"mmseqs" if auto_cluster and not cluster else None}))
        st.caption("CataPro 下载命令和校验值见升级说明。SABIO-RK 保留原始 SBML，并按序列、底物和参数逐项映射。")
    with st.expander("提取特征与训练"):
        dataset=st.text_input("数据集目录",key="train_dataset")
        device=st.selectbox("特征计算设备",["cuda","cpu"],key="feature_device")
        if st.button("提取冻结 ESM-2 特征"):_do(lambda:_submit("features",{"dataset":dataset,"device":device}))
        features=st.text_input("已完成的特征目录",key="train_features")
        conditions=st.checkbox("使用温度和 pH（缺失条件的样本不进入本次训练）")
        if st.button("训练并评估 Ridge / ExtraTrees"):
            _do(lambda:_submit("train",{"dataset":dataset,"features":features,"conditions":conditions}))
    with st.expander("登记模型版本"):
        finished=[r for r in Store().tasks(context()["project_id"]) if r["state"]=="completed" and json.loads(r["spec"])["kind"]=="train"]
        if finished:
            trained=st.selectbox("已完成的训练任务",[r["id"] for r in finished])
            if st.button("登记此训练模型"):
                result=_do(lambda:register_task_model(trained,context()["project_id"]))
                if result:st.success(str(result))
        directory=st.text_input("模型目录",key="model_register")
        if st.button("验证并登记"):
            result=_do(lambda:register_model(directory,context()["project_id"]))
            if result:st.success(str(result))
    if models:
        with st.expander("用自主模型预测 kcat / Km",expanded=True):
            from ..property_screen import parse_substrates
            choices={r["card"]["model_id"]:r for r in models}
            selected=choices[st.selectbox("模型版本",list(choices),key="kinetic_model")]
            enzymes=Store().records('enzyme',context()['project_id'])
            if enzymes:
                with st.expander('从项目中选择酶与构建'):
                    labels={r['enzyme_id']:r['name']+' · '+str(len(r['sequence']))+' aa · '+r['construct_id'] for r in enzymes}
                    chosen=st.multiselect('待预测的酶记录',list(labels),format_func=labels.__getitem__,key='kinetic_enzyme_records')
                    if st.button('载入所选序列',disabled=not chosen):
                        st.session_state['kinetic_fasta']=''.join('>'+r['enzyme_id']+'\n'+r['sequence']+'\n' for r in enzymes if r['enzyme_id'] in chosen)
            sequence=st.text_area("待预测序列 FASTA",key="kinetic_fasta")
            chemicals=st.text_area("底物名称=SMILES",value="octanoic_acid=CCCCCCCC(=O)O",key="kinetic_smiles")
            conditional=selected['card']['conditions']
            if conditional:
                temperature=st.number_input("测定温度（K）",min_value=200.,max_value=500.,value=313.15,key="kinetic_temperature")
                ph=st.number_input("测定 pH",min_value=0.,max_value=14.,value=7.,key="kinetic_ph")
            else:
                temperature=ph=None
                st.caption('此版本为无条件模型，温度和 pH 不参与预测。')
            if st.button("提交动力学预测"):
                def infer():
                    executor=selected.get("executor","local")
                    identifier=submit("kinetics",{"model_dir":selected["path"],"text":sequence,
                        "substrates":parse_substrates(chemicals.splitlines()),"temperature_K":temperature,"pH":ph,
                        "device":"cuda" if executor!="local" else "cpu"},
                        project_id=context()["project_id"],family_id=context()["family_id"],executor=executor)
                    st.success("预测任务已保存："+identifier)
                    return identifier
                _do(infer)
        for model in models:
            card=model["card"]
            with st.expander(card['endpoint']+' 评估详情 · '+card['model_id']):
                st.json({"指标":card["metrics"],"改善区间":card["promotion"],"校准":card["calibration"]})
    else:st.info("尚未登记自主训练模型。已有 SWI、TmProt、EnzRank 仍在性质预测页独立使用。")
    reports=kinetic_reports(project_id=context()['project_id'])
    if reports:
        with st.expander('动力学预测结果',expanded=True):
            report=reports[st.selectbox('已完成的动力学预测',list(reports),key='kinetic_report')]
            payload=json.loads((report/'predictions.json').read_text(encoding='utf-8'))
            rows=[]
            for item in payload['predictions']:
                rows.append({'候选':', '.join(item.get('candidate_ids',[item['enzyme_id']])),
                             '底物':', '.join(item.get('substrate_names',[item['reaction_id']])),
                             '指标':item['endpoint'],'预测值':item['value'],'单位':item['unit'],
                             '适用范围':{'predicted':'训练域内','out_of_domain':'训练域尚未确认'}.get(item['status'],item['status']),
                             '模型版本':item['model_id']})
            st.dataframe(pd.DataFrame(rows),hide_index=True,width='stretch')
            st.caption('结果是动力学预测；不能换算为酯合成收率或合成/水解比。训练域尚未确认的序列不发布校准区间。')
            for filename,mime in [('predictions.csv','text/csv'),('predictions.json','application/json')]:
                st.download_button('下载 '+filename,(report/filename).read_bytes(),filename,mime,key='kinetic_download_'+filename)

def mechanism_page():
    from .mechanism_plan import MechanismBatch,freeze_batch
    st.subheader("机制计算")
    st.caption("首批为游离酶：亲本加 3 条候选，每个 3 次独立重复，每次生产模拟 10 ns，总计 120 ns；平衡阶段另计。")
    st.write("流程：构建确认 → 装配与参数检查 → 最小化和平衡 → 分段 MD → 轨迹质量与通道、水/醇分析。")
    preparations=Store().records('preparation',context()['project_id'])
    if preparations:
        prepared=preparations[-1]
        st.subheader('当前建模与参数检查')
        if prepared.get('structures'):
            rows=[]
            for row in prepared['structures']:
                state={'queued':'排队中','running':'预测中','failed':'失败','completed':'已完成'}.get(row['state'],row['state'])
                rows.append({'构建':row['name'],'任务':state,'结构复核':('通过' if row.get('structure_review_pass') else '待通过') if row['state']=='completed' else '待计算',
                             '平均 pLDDT':row.get('plddt'),'N 端 pLDDT':row.get('terminal_plddt'),'Ser–His 距离（Å）':row.get('ser_his_A')})
            st.dataframe(pd.DataFrame(rows),hide_index=True,width='stretch')
            st.caption('结构复核采用完整序列、重原子、原子冲突、催化几何与置信度检查；正式模拟还需装配和溶液平衡通过。')
        if prepared.get('parameters'):
            with st.expander(f"配体与缓冲组分：{len(prepared['parameters'])} 种状态",expanded=True):
                st.dataframe(pd.DataFrame([{'分子状态':r['name'],'形式电荷':r['formal_charge'],'部分电荷总和':r['partial_charge_sum'],
                     '原子数':r['atoms'],'参数检查':'数值检查通过' if r['status']=='numerical_checks_passed' else '待通过'} for r in prepared['parameters']]),hide_index=True,width='stretch')
                st.caption('GAFF2 / AM1-BCC 参数数值检查与水溶液验证分别记录。柠檬酸位置异构微状态的比例仍属于模拟假设。')
        if prepared.get('solutions'):
            st.dataframe(pd.DataFrame([{'溶液对照':r['name'],'任务':{'queued':'排队中','running':'计算中','failed':'失败','completed':'已完成'}.get(r['state'],r['state']),
                '密度（g/mL）':r.get('density_g_ml'),'密度分段漂移':r.get('density_drift')} for r in prepared['solutions']]),hide_index=True,width='stretch')
            st.caption('上述为不含蛋白的溶液测试，时长单独记录；120 ns 正式酶模拟的完成量另行统计。')
        if prepared.get('assemblies'):
            st.write('二聚体装配复核')
            st.dataframe(pd.DataFrame([{'构建':r['name'],'装配检查':{'passed':'通过','blocked':'待修正'}.get(r.get('status'),{'running':'计算中','queued':'排队中','failed':'失败'}.get(r['state'],r['state'])),
                '链间严重碰撞数':r.get('clashes'),'核心骨架位移（Å）':r.get('core_rms_displacement_A')} for r in prepared['assemblies']]),hide_index=True,width='stretch')
            st.caption('按公开生物学装配放置两条链，固定核心骨架约束后修正末端碰撞；这一步不等同于水溶液平衡。')
        if prepared.get('bulk_analysis'):
            with st.expander('溶液分布与接触分析'):
                st.dataframe(pd.DataFrame([{'溶液':r['name'],'检查':r.get('status',r['state']),
                    '柠檬酸分子数':r.get('citrate_count'),'平均最大接触簇':r.get('citrate_largest_cluster_mean'),
                    '观察到的最大接触簇':r.get('citrate_largest_cluster_max')} for r in prepared['bulk_analysis']]),hide_index=True,width='stretch')
                st.caption('接触按周期边界下的原子距离计算；短轨迹只报告观察到的分布，不据此判定聚集平衡或结合常数。')
        if prepared.get('parameter_bundle') and Path(prepared['parameter_bundle']).is_file():
            st.download_button('下载已检查的分子参数',Path(prepared['parameter_bundle']).read_bytes(),'parameters.zip','application/zip')
        st.info('正式模拟仍需逐构建完成蛋白质子化设定、显式溶液装配和分阶段平衡检查。')
    for record in Store().records("readiness",context()["project_id"]):
        st.write(f"完整构建已映射：{record.get('mapped_designs',0)} 条；正式模拟候选已放行：{record.get('selected_production_candidates',0)} 条。")
        reasons={"full_construct_coordinates_include_residues_1_11":"补全并检查 N 端 1–11 位残基",
                 "independent_candidate_structure_review":"候选完整构建的独立结构复核",
                 "validated_citrate_and_ligand_parameters":"柠檬酸及配体参数验证",
                 "geometry_and_equilibration_acceptance":"催化几何与平衡检查通过"}
        if preparations:
            with st.expander('早期准备记录（保留原始缺项）'):
                for item in record.get('blocking_items',[]):st.write(reasons.get(item,item))
        else:
            for item in record.get("blocking_items",[]):st.warning(reasons.get(item,item))
    from .artifacts import result_roots
    bundles={}
    for root in result_roots(context()['project_id']):
        for path in root.glob('mechanism/parameters.zip'):bundles[str(path)]=path
    if bundles:
        with st.expander('已完成的参数任务下载'):
            keys=list(bundles)
            selected=st.selectbox('参数任务',keys,format_func=lambda key:next((part for part in Path(key).parts if part.startswith('task_')),'分子参数'))
            st.download_button('下载参数文件与检查报告',bundles[selected].read_bytes(),'parameters.zip','application/zip',key='parameter_task_download')
    with st.expander("冻结批次与计算预算",expanded=True):
        manifest=st.file_uploader("候选与反应批次 JSON",type=["json"],key="mechanism_batch")
        ns=st.number_input("每次生产模拟长度（ns）",min_value=0.1,max_value=1000.,value=10.,step=1.)
        replicas=st.number_input("每个序列的独立重复",1,20,3)
        st.metric("亲本加 3 条候选的总生产时长",f"{4*replicas*ns:g} ns")
        if manifest and st.button("冻结本批次"):
            def save():
                raw=json.loads(manifest.getvalue())
                raw.update(production_ns=ns,replicas=replicas,approved_budget_ns=len(raw["candidates"])*ns*replicas,seeds=[])
                batch=MechanismBatch.model_validate(raw)
                path=state_root()/"batches"/(batch.batch_id+".json")
                result=freeze_batch(batch.model_dump(),path)
                Store().record("mechanism_batch",batch.batch_id,context()["project_id"],batch.model_dump())
                return result
            result=_do(save)
            if result:st.success(str(result))
    with st.expander('生成分子参数'):
        st.caption('每行输入名称、SMILES 和形式电荷，以空格分隔。适用于已明确质子化状态的 H/C/N/O 小分子。')
        molecule_text=st.text_area('分子清单',value='ethanol CCO 0\noctanoic_acid CCCCCCCC(=O)O 0\noctanoate CCCCCCCC(=O)[O-] -1',key='parameter_molecules')
        if st.button('生成 GAFF2 参数并检查'):
            def submit_parameters():
                molecules=[]
                for line in molecule_text.splitlines():
                    if not line.strip():continue
                    fields=line.split()
                    if len(fields)!=3:raise ValueError('每行需要名称、SMILES 和形式电荷三项')
                    molecules.append({'name':fields[0],'smiles':fields[1],'formal_charge':int(fields[2])})
                request={'version':'1.0.0','forcefield':'GAFF2','charge_method':'AM1-BCC','molecules':molecules}
                path=state_root()/'ui_inputs'/(uid('molecules')+'.json');path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(json.dumps(request),encoding='utf-8')
                return _submit('mechanism',{'action':'parameters','input':str(path)})
            _do(submit_parameters)
    with st.expander("准备与运行"):
        st.caption("输入已确认的构建和参数清单；缺失残基、未支持的金属/共价状态或未通过的质量检查会列出具体缺项。")
        action=st.selectbox("执行阶段",["parameters","solution","solution-analyze","relax-assembly","check","prepare","run","batch","analyze"],
                            format_func=lambda x:{"parameters":"生成配体参数与数值检查","solution":"溶液密度与浓度检查","solution-analyze":"溶液分布与接触分析","relax-assembly":"修正装配碰撞并复核","check":"检查准备状态","prepare":"装配与平衡","run":"分段生产模拟","batch":"执行已冻结的独立重复批次","analyze":"通道和水/醇分析"}[x])
        upload=st.file_uploader("此阶段输入 JSON",type=["json"],key="mechanism_input")
        if st.button("提交机制任务",disabled=upload is None):
            _do(lambda:_submit("mechanism",{"action":action,**json.loads(upload.getvalue())}))
    st.info("120 ns 批次用于流程验证与初步采样。MSM/TPT 速率需独立通过转移事件、隐含时间尺度和 Chapman–Kolmogorov 检查。")
