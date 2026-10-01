-- Proteins with an AMR annotation.
SELECT DISTINCT a.id
FROM amr a
JOIN proteins p ON p.id = a.id
