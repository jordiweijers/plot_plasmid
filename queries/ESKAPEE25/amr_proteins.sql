-- Proteins on the contig with an AMR annotation. Runs once per contig.
SELECT DISTINCT a.id
FROM amr a
JOIN proteins p ON p.id = a.id
WHERE p.contig_id = @contig_id
