"""Project-scoped mechanism results in the existing scientific workbench."""
from __future__ import annotations
import json
from datetime import datetime
from pathlib import Path
import pandas as pd
import streamlit as st
from .store import Store
from .datasets import file_hash
from .mechanism_results import load_result

SOLVENTS={'water':'水','ethanol':'乙醇','octanoic_acid':'辛酸'}

@st.cache_data(show_spinner=False)
def cached_result(path,modified,size,dependencies,quality_record):
    return load_result(path,quality_record)


def render_results(project_id):
    store=Store();choices={};descriptions={};reviews={r['task_id']:r for r in store.records('mechanism_quality',project_id)}
    for task in store.tasks(project_id):
        spec=json.loads(task['spec'])
        if spec['kind']!='mechanism' or spec['inputs'].get('action')!='batch':continue
        result=json.loads(task.get('result') or '{}')
        root=result.get('local_output')
        if not root:root=str(store.path.parent/'tasks'/task['id']/'output')
        path=Path(root)/'mechanism/batch_report.json'
        if path.is_file():
            choices[task['id']]=path
            metadata=result.get('remote') or result
            date=datetime.fromisoformat(task['created']).astimezone().strftime('%m-%d %H:%M')
            descriptions[task['id']]=f"{date} · {metadata.get('production_ns','—')} ns · {metadata.get('replicas_completed','—')} 条轨迹"
    if not choices:return False
    st.subheader('生产模拟结果')
    selected=st.selectbox('查看已完成批次',list(choices),format_func=descriptions.get,key='mechanism_result_batch')
    path=choices[selected];quality_record=reviews.get(selected)
    try:
        stat=path.stat()
        # Recheck even when the parsed report is cached.
        if quality_record and file_hash(quality_record['report_path'])!=quality_record['sha256']:
            raise ValueError('Quality report checksum mismatch')
        dependencies=tuple((str(p),file_hash(p)) for p in sorted([path.parent/'batch_identity.json',
            *path.parent.glob('*/preflight.json'),*path.parent.glob('*/trajectory_report.json')]) if p.is_file())
        data=cached_result(str(path),stat.st_mtime_ns,stat.st_size,dependencies,quality_record)
    except (OSError,ValueError,KeyError,TypeError) as exc:
        st.error('结果文件未通过读取检查，请重新同步该任务：'+str(exc));return True
    columns=st.columns(4)
    columns[0].metric('生产模拟',f"{data['production_ns']:g} ns")
    columns[1].metric('完成轨迹',str(data['replicas']))
    columns[2].metric('位点分析完整',f"{data['analyzed_replicas']}/{data['replicas']}")
    columns[3].metric('轨迹数值质检通过',f"{data['quality_passed']}/{data['replicas']}")
    st.caption('这里展示模拟描述量。先在每条轨迹内平均两个或多个位点，再比较独立重复；亚基不另算独立样本。')
    labels=(quality_record or {}).get('candidate_labels',{})
    def display(frame):
        frame=frame.copy()
        if 'candidate_id' in frame:frame['candidate_id']=frame['candidate_id'].map(lambda x:labels.get(x,x))
        if 'solvent' in frame:frame['solvent']=frame['solvent'].map(lambda x:SOLVENTS.get(x,x))
        return frame
    compare,sites,qc,files=st.tabs(['候选与重复比较','逐位点结果','轨迹质检','报告下载'])
    with compare:
        frame=pd.DataFrame(data['candidates']);replicates=pd.DataFrame(data['replica_values'])
        if frame.empty:st.info('尚无位点齐全的轨迹分析；可在逐位点结果中查看已有数据。')
        else:
            solvent=st.selectbox('比较的分子',list(dict.fromkeys(frame.solvent)),format_func=lambda x:SOLVENTS.get(x,x),key='mechanism_result_solvent')
            subset=display(frame[frame.solvent==solvent])
            columns={'candidate_id':'构建','replicates':'独立重复数',
                'mean_molecules_mean':'位点平均分子数','mean_molecules_min':'重复最小值','mean_molecules_max':'重复最大值',
                'occupied_fraction_mean':'有占据的帧比例'}
            st.dataframe(subset[list(columns)].rename(columns=columns),hide_index=True,width='stretch',
                column_config={'构建':st.column_config.TextColumn(width='medium')})
            st.caption(f"位点半径：{subset['radius_A'].iloc[0]:g} Å；进入次数可在下方逐重复指标和逐位点结果中查看。")
            import altair as alt
            points=display(replicates[replicates.solvent==solvent])
            metric=st.selectbox('逐重复指标',['mean_molecules','occupied_fraction','entry_events'],
                format_func=lambda x:{'mean_molecules':'位点平均分子数','occupied_fraction':'有占据的帧比例','entry_events':'每位点进入次数'}[x],key='mechanism_plot_metric')
            chart=alt.Chart(points).mark_circle(size=95,opacity=.8).encode(
                x=alt.X('candidate_id:N',title='构建',sort=None,axis=alt.Axis(labelAngle=0)),
                y=alt.Y(metric+':Q',title={'mean_molecules':'平均分子数','occupied_fraction':'占据帧比例','entry_events':'观察到的进入次数'}[metric],scale=alt.Scale(zero=True)),
                color=alt.Color('candidate_id:N',legend=None,scale=alt.Scale(range=['#12665D','#3E6C99','#B58538','#876788'])),
                tooltip=[alt.Tooltip('candidate_id:N',title='构建'),alt.Tooltip('replica_id:N',title='独立重复'),alt.Tooltip(metric+':Q',format='.3f')]).properties(height=260)
            st.altair_chart(chart,width='stretch')
            st.caption('每个点是一条独立重复轨迹；表中范围是重复间的最小值和最大值，未换算为酶活或合成/水解比。')
    with sites:
        frame=display(pd.DataFrame(data['sites']))
        columns={'candidate_id':'构建','replica_id':'重复','site':'位点','solvent':'分子','frames':'帧数','mean_molecules':'平均分子数',
            'occupied_fraction':'占据帧比例','entry_events':'进入次数','complete_residence_events':'完整驻留事件数',
            'complete_event_observed_median_ps':'完整事件观测中位时长（ps）','left_censored_events':'轨迹开始前已进入的事件',
            'right_censored_events':'轨迹结束时尚未离开的事件','channel_sample_frames':'通道采样帧数',
            'path_found_fraction':'找到路径的采样帧比例','path_bottleneck_mean_A':'已找到路径的平均瓶颈（Å）'}
        if not frame.empty:st.dataframe(frame[list(columns)].rename(columns=columns),hide_index=True,width='stretch')
        st.caption('轨迹边界截断的驻留事件单独计数。未找到通道的帧保留为未找到，不按零瓶颈参与均值；通道结果仍需网格与探针敏感性检查。')
    with qc:
        frame=display(pd.DataFrame(data['quality']))
        frame['quality_status']=frame.quality_status.map({'passed':'数值检查通过','requires_review':'需复核','unreviewed':'尚未质检'})
        columns={'candidate_id':'构建','replica_id':'重复','production_ns':'生产时长（ns）','analysis_complete':'位点分析完整',
            'quality_status':'质检状态','core_rmsd_mean_A':'核心平均 RMSD（Å）','core_rmsd_max_A':'核心最大 RMSD（Å）','issues':'发现的问题'}
        st.dataframe(frame[[c for c in columns if c in frame]].rename(columns=columns),hide_index=True,width='stretch')
        st.caption('数值检查、结构变化和机制采样充分性分别判断。RMSD 为相对各构建平衡结构的核心 Cα 位移；生产阶段催化距离作为描述量保留。')
        sampling=data['sampling_status']
        st.write('MSM/TPT 采样评估：'+('初步采样，尚未通过收敛检验' if sampling=='pilot; MSM/TPT convergence not established' else sampling))
        if data['assumptions']:
            with st.expander('本批次的质子化与采样假设'):
                for assumption in data['assumptions']:st.write(assumption)
    with files:
        for key,label in [('candidates','候选比较'),('replica_values','逐重复结果'),('sites','逐位点结果'),('quality','轨迹质检')]:
            frame=display(pd.DataFrame(data[key]));blob=frame.to_csv(index=False).encode('utf-8-sig')
            st.download_button('下载'+label+' CSV',blob,file_name=key+'.csv',mime='text/csv',key='mechanism_download_'+key)
        st.download_button('下载原始批次分析报告',path.read_bytes(),'batch_report.json','application/json',key='mechanism_report_json')
        if data['quality_report_path']:
            st.download_button('下载完整结构质检报告',Path(data['quality_report_path']).read_bytes(),'structural_review.json','application/json',key='mechanism_quality_json')
    return True


def render_structure_reviews(project_id):
    records=Store().records('structure_batch_review',project_id)
    if not records:return
    record=records[-1]
    try:
        path=Path(record['report_path'])
        if file_hash(path)!=record['sha256']:raise ValueError('Structure-review checksum mismatch')
        report=json.loads(path.read_text(encoding='utf-8'))
        rows=[{'构建':r['candidate_id'],'结构复核':'通过' if r.get('structure_review_pass') else '需复核',
               '平均 pLDDT':r.get('mean_plddt'),'N 端平均 pLDDT':r.get('terminal_plddt'),
               '问题':'; '.join(r.get('blocking_issues',[]))} for r in report['rows']]
    except (OSError,ValueError,KeyError,TypeError) as exc:
        st.error('完整构建复核报告读取失败：'+str(exc));return
    with st.expander(f"完整构建结构复核：{sum(r.get('structure_review_pass',False) for r in report['rows'])}/{report['total']} 通过"):
        st.dataframe(pd.DataFrame(rows),hide_index=True,width='stretch')
        st.caption('完整序列、重原子、碰撞、催化位点和置信度分别检查；该结果不替代蛋白装配、平衡或实验验证。')
        st.download_button('下载完整构建复核报告',path.read_bytes(),'structure_review.json','application/json',key='full_construct_review_download')
