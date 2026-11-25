"""Exercise each original workbench page, not only the initially selected tab."""
from streamlit.testing.v1 import AppTest
from bur01_structmine.config import project_root


def test_all_workbench_pages_render_in_chinese(tmp_path,monkeypatch):
    monkeypatch.setenv('AQUA_TEST_MODE','1')
    for name in ['RUN','PREDICTION','NR_RUN','DESIGN','PROPERTY']:
        monkeypatch.setenv('AQUA_'+name+'_ROOT',str(tmp_path/name.lower()))
    monkeypatch.setenv('AQUA_REMOTE_ACCESS_PATH',str(tmp_path/'remote_access.json'))
    app=AppTest.from_file(str(project_root()/'ui/app.py'),default_timeout=30).run()
    tabs=[tab.label for tab in app.tabs]
    assert len(tabs)==13
    for label in tabs:
        app.session_state['structmine_tab']=label
        app.run()
        assert not app.exception,(label,[e.message for e in app.exception])
        for element in list(app.markdown)+list(app.caption)+list(app.error):
            assert '\ufffd' not in element.value,(label,element.value)
