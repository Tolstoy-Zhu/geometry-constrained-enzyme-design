"""Show engineering acceptance independently of research validation."""
import json
from datetime import datetime
from pathlib import Path
import pandas as pd
import streamlit as st
from .datasets import file_hash
from .store import Store


def render_acceptance(project_id):
    records=Store().records('platform_acceptance',project_id)
    if not records:return
    record=records[-1]
    with st.expander('平台整体验收',expanded=True):
        try:
            path=Path(record['report_path'])
            if file_hash(path)!=record['sha256']:raise ValueError('验收报告校验值不一致')
            report=json.loads(path.read_text(encoding='utf-8'))
            if report['project_id']!=project_id:raise ValueError('报告不属于当前项目')
        except (OSError,ValueError,KeyError) as exc:st.error(str(exc));return
        passed=sum(r['status']=='passed' for r in report['checks'])
        st.write(f"工程检查：{passed}/{len(report['checks'])} 项通过")
        stamp=datetime.fromisoformat(report['created_at']).astimezone().strftime('%Y-%m-%d %H:%M')
        st.caption('验收时间：'+stamp+'。此报告记录该时点的状态；后续代码、数据或环境变化需重新验收。')
        summaries={'families':lambda d:'Bur01、CALB 和糖酶分别使用家族规则',
            'records':lambda d:f"{d['enzymes']} 条酶/构建，{d['reactions']} 组反应，{d['assays']} 条测定记录",
            'structures':lambda d:f"{d['passed']}/{d['reviewed']} 个完整构建通过结构复核",
            'properties':lambda d:f"{d['readable_reports']} 份性质预测报告可读取",
            'kinetics':lambda d:f"{d['readable_reports']} 份动力学预测报告可读取",
            'models':lambda d:'、'.join(r['endpoint'] for r in d)+' 的模型版本和分项评估已登记',
            'production':lambda d:'；'.join(f"{r['replicas']} 条轨迹、{r['production_ns']:g} ns，逐位点分析和数值质检齐全" for r in d),
            'mechanism_assessment':lambda d:'；'.join(f"{r['parameter_evaluations']} 次参数检查、{r['sampling_groups']} 组采样诊断" for r in d),
            'history':lambda d:f"{d['unchanged_history_manifests']} 份历史记录保持原样",
            'feedback':lambda d:'重复导入不重号，冲突不覆盖；模拟值和未测值均被训练输入排除',
            'web':lambda d:f"{d['pages']} 个页面已检查",
            'remote_recovery':lambda d:f"原任务编号恢复，{d['unchanged_saved_frames']} 份已完成帧未被重写"}
        table=[]
        for row in report['checks']:
            detail=summaries[row['id']](row['detail']) if row['status']=='passed' and row['id'] in summaries else str(row['detail'])
            if detail=='Parameter and sampling assessment is not complete':detail='最新一轮参数与采样检查尚未完成'
            table.append({'检查项':row['title'],'状态':'通过' if row['status']=='passed' else '需处理','结果':detail})
        st.dataframe(pd.DataFrame(table),hide_index=True,width='stretch')
        research=report['research_validation']
        st.write(f"本项目精确实测记录：{research['experimental_assays']} 条；机制速率尚未建立；完整 NR 验收单独保留。")
        st.caption('实验回填测试使用独立软件测试索引；模拟值和未测记录均被实测训练集排除。')
        st.download_button('下载平台验收报告',path.read_bytes(),'platform_acceptance.json','application/json',key='platform_acceptance_download')
