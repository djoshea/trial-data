function [td, tdi, R] = loadTrialsIntoTrialData(dateStr, subject, protocol)
    pathMgr = FilePathManager();
    pathMgr.protocol = protocol;
    pathMgr.subject = subject;
    pathMgr.date = datenum(dateStr);

    tdl = TrialDataLoader(pathMgr);

    R = tdl.loadAllTrialsIntoStructArray();
    tdi = MatUdpTrialDataInterface(R);
    td = TrialDataConditionAlign(tdi);
end
