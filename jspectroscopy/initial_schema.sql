create table fits (
    detnum INTEGER NOT NULL,
    x    INTEGER NOT NULL,
    y    INTEGER NOT NULL,
    group_num INTEGER NOT NULL,
    success INTEGER NOT NULL,
    params_result BLOB,
    params_result_errors BLOB,
    cov BLOB,
    fit_bounds BLOB,
    redchisqr REAL NOT NULL,
    PRIMARY KEY( detnum, x, y, group_num )
);

-- create table fits (
--     x    INTEGER NOT NULL,
--     y    INTEGER NOT NULL,
--     fit_id INTEGER NOT NULL,
--     successful_fit INTEGER NOT NULL,
--     npeaks INTEGER NOT NULL,
--     last_attempt INTEGER NOT NULL,
--     params_guess BLOB,
--     fit_bounds BLOB,
--     peak_guesses BLOB,
--     model BLOB,
--     PRIMARY KEY( x, y, fit_id )
-- );



create table metadata (
    dets_used STRING NOT NULL,
    xdim INTEGER NOT NULL,
    ydim INTEGER NOT NULL,
    peak_types BLOB,
    constrain_det_params BLOB,
    timestamp TEXT,
    name TEXT
);
