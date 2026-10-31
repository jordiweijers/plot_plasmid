-- Proteins with an AMR annotation in kes2021.
-- A category query runs on the whole database and returns one column of IDs; here protein IDs, so use it with --id-type protein.
SELECT DISTINCT a.id
FROM amr a
JOIN proteins p ON p.id = a.id
