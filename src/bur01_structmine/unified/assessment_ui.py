"""Parameter sensitivity and sampling checks alongside the original MD results."""
from __future__ import annotations
import json
import pandas as pd
import streamlit as st
from .assessment_results import assessment_tasks,load_assessment,submit_assessment

STATUSES={'grid_sensitive':'对网格设置敏感','edge_clearance_requires_review':'点间路径存在碰撞，需复核',
    'no_paths_on_sampled_frames':'采样帧中未找到路径','grid_agreement_on_sampled_frames':'网格差异在设定阈值内'}
REASONS={'too_few_independent_replicas':'独立重复数不足','too_few_bidirectional_occupancy_transitions':'双向占据转移事件不足',
    'insufficient_empty_or_occupied_state_frames':'空位或占据状态的帧数不足',
    'no_resolved_implied_timescale_plateau':'未见可分辨的时间尺度平台',
    'held_out_chapman_kolmogorov_not_passed':'留一重复的 CK 检查未通过',
    'occupancy_states_not_validated_for_transport_TPT':'占据状态尚未验证为分子输运状态'}
SOLVENTS={'water':'水','ethanol':'乙醇','octanoic_acid':'辛酸'}


def render_assessment(project_id,production_task_id,production_report,labels,channel_tab,sampling_tab):
    rows=assessment_tasks(project_id,production_task_id)
    completed=next((r for r in rows if r['state']=='completed'),None)
    if not completed:
        with channel_tab:
            if rows:
                state=rows[0]['state'];st.info('参数与采样检查：'+{'running':'正在分析现有轨迹','queued':'排队中','cancelled':'已停止，可在项目与任务恢复','failed':'分析失败，请查看任务日志','blocked':'等待重新连接确认'}.get(state,state))
            else:
                st.info('此批次尚未运行通道参数与采样评估。')
                if st.button('检查现有轨迹的通道参数与采样',key='start_mechanism_assessment'):
                    try:
                        submit_assessment(project_id,production_task_id);st.rerun()
                    except (ValueError,OSError,RuntimeError,KeyError) as exc:st.error(str(exc))
            st.caption('仅分析已有轨迹，默认比较 3 档网格、4 档探针和 2 个网格起点；不会追加分子动力学时长。')
        with sampling_tab:st.info('采样评估完成并同步后，在此显示逐构建、逐位点的诊断结果。')
        return
    try:report,root=load_assessment(completed,production_task_id,production_report)
    except (ValueError,OSError,KeyError,TypeError) as exc:
        with channel_tab:st.error('评估报告读取检查失败：'+str(exc))
        return
    def label(candidate):return labels.get(candidate,candidate)
    with channel_tab:
        protocol=report['protocol']
        if rows and rows[0]['state'] in {'queued','running'}:
            st.info('更新后的检查正在运行；下方暂时显示上一份已完成报告。')
        if protocol.get('channel_method','point_grid_v1')=='point_grid_v1':
            st.warning('这是原网格点算法的诊断报告。逐线段检查发现的碰撞路径不视为通过；修正版结果单独保存。')
            if st.button('使用修正版重新检查已有轨迹',disabled=bool(rows and rows[0]['state'] in {'queued','running'}),key='upgrade_mechanism_assessment'):
                try:
                    submit_assessment(project_id,production_task_id);st.rerun()
                except (ValueError,OSError,RuntimeError,KeyError) as exc:st.error(str(exc))
        else:
            st.caption('当前版本：连续线段间隙已计入路径搜索，并对返回路径再次独立复核。')
        st.caption('网格间距：'+', '.join(map(str,protocol['spacings_A']))+' Å；探针半径：'+', '.join(map(str,protocol['probes_A']))+' Å。另比较网格起点偏移，并验证网格点之间整段路径的原子间隙。')
        frame=pd.DataFrame(report['channel']['rows'])
        probe=st.selectbox('探针半径（Å）',sorted(frame.probe_radius_A.unique()),key='assess_probe')
        site=st.selectbox('通道位点',list(dict.fromkeys(frame.site)),key='assess_channel_site')
        subset=frame[(frame.probe_radius_A==probe)&(frame.site==site)].copy()
        subset['candidate_id']=subset.candidate_id.map(label);subset['replica_id']=subset.replica_id.str.rsplit('_',n=1).str[-1]
        subset['status']=subset.status.map(STATUSES)
        columns={'candidate_id':'构建','replica_id':'重复','sampled_frames':'采样帧',
            'grid_disagreement_fraction':'网格判定不一致比例','matched_bottleneck_range_max_A':'共同有路径帧的最大瓶颈差（Å）',
            'edge_clearance_failures':'点间碰撞次数','status':'检查结果'}
        st.dataframe(subset[list(columns)].rename(columns=columns),hide_index=True,width='stretch')
        st.write(f"共完成 {report['channel']['frame_parameter_evaluations']:,} 次位点与参数组合检查。")
        st.caption('不一致比例按同一轨迹的采样帧计算；无共同有路径帧时不评估瓶颈差。网格设置和两个亚基不增加独立重复数。探针半径尚未标定为底物有效尺寸；未找到路径不等于证实通道关闭。')
        st.download_button('下载通道参数检查 CSV',(root/'channel_parameter_summary.csv').read_bytes(),'channel_parameter_summary.csv','text/csv',key='assessment_channel_csv')
        st.download_button('下载全部采样帧与参数 CSV',(root/'channel_parameter_frames.csv').read_bytes(),'channel_parameter_frames.csv','text/csv',key='assessment_frames_csv')
    with sampling_tab:
        groups=report['sampling']['groups']
        solvent=st.selectbox('评估分子',list(dict.fromkeys(r['solvent'] for r in groups)),format_func=lambda x:SOLVENTS.get(x,x),key='assess_solvent')
        site=st.selectbox('采样位点',list(dict.fromkeys(r['site'] for r in groups)),key='assess_sampling_site')
        selected=[r for r in groups if r['solvent']==solvent and r['site']==site]
        table=[{'构建':label(r['candidate_id']),'独立重复':r['independent_replicas'],
            '双向转移最少次数':min(min(v['empty_to_occupied'],v['occupied_to_empty']) for v in r['replicas']),
            '时间尺度平台':'通过' if r['its_pass'] else '未通过','留一重复 CK':'通过' if r['ck_pass'] else '未通过',
            '尚缺条件':'；'.join(REASONS.get(x,x) for x in r['reasons'])} for r in selected]
        st.dataframe(pd.DataFrame(table),hide_index=True,width='stretch')
        st.info('当前状态：尚未建立可靠的 MSM/TPT 机制速率。这里诊断的是位点“空位/有分子占据”两种状态，不能直接换算成单分子的进入速率、酶活或合成/水解比。')
        candidate=st.selectbox('查看采样诊断详情',[r['candidate_id'] for r in selected],format_func=label,key='assess_candidate')
        detail=next(r for r in selected if r['candidate_id']==candidate)
        replicas=[{'重复':r['replica_id'].rsplit('_',1)[-1],'空位帧':r['empty_frames'],'占据帧':r['occupied_frames'],
                   '空位→占据':r['empty_to_occupied'],'占据→空位':r['occupied_to_empty'],
                   '前半段平均分子数':r['first_half_mean'],'后半段平均分子数':r['second_half_mean']} for r in detail['replicas']]
        st.dataframe(pd.DataFrame(replicas),hide_index=True,width='stretch')
        its=pd.DataFrame([{'滞后时间（ps）':r['lag_ps'],'隐含时间尺度（ps）':r['timescale_ps'],
            '状态':{'estimated':'已估计','unobserved_state':'存在未采到的状态','disconnected':'未连通','unresolved_or_oscillatory':'不可分辨或振荡'}.get(r['status'],r['status'])} for r in detail['implied_timescales']])
        st.dataframe(its,hide_index=True,width='stretch')
        with st.expander('留一重复的 Chapman–Kolmogorov 检查'):
            ck=[{'留出的重复':r['held_out_replica'].rsplit('_',1)[-1],'滞后时间（ps）':r['lag_ps'],
                '最大概率差':r['max_absolute_error'],'状态':'已计算' if r['status']=='evaluated' else '状态计数不足'} for r in detail['held_out_ck']]
            st.dataframe(pd.DataFrame(ck),hide_index=True,width='stretch')
        st.caption('三个随机速度重复共享各构建的平衡起始结构。检查阈值在本次分析前固定，是初步筛查规则；通过占据诊断也不代表机制状态或平衡采样已验证。')
        st.download_button('下载采样评估 CSV',(root/'sampling_summary.csv').read_bytes(),'sampling_summary.csv','text/csv',key='assessment_sampling_csv')
        st.download_button('下载完整参数与采样评估',(root/'assessment_report.json').read_bytes(),'assessment_report.json','application/json',key='assessment_json')
