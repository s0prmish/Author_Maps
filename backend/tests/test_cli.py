"""Tests for the CLI."""
import os
import click
import pytest
from click.testing import CliRunner
from authormaps import cli

class TestCli:

    def test_get_publications(self):
        runner = CliRunner()

        result = runner.invoke(cli.get_publications, "'Birkenbihl Colin'")
        assert result
        assert result.exit_code == 0
        assert "33344750" in result.output
        assert "33285634" in result.output
        assert "32843907" in result.output
        assert "31993440" in result.output

        result = runner.invoke(cli.get_publications, "'Yalchyk'")
        assert result
        assert result.exit_code == 0
        assert result.output == "No publications for 'Yalchyk'\n"  # which is sad

    def test_get_coauthors(self):
        runner = CliRunner()

        result = runner.invoke(cli.get_coauthors, "'Birkenbihl Colin'")
        assert result
        assert result.exit_code == 0
        assert "Domingo-Fernandez Daniel" in result.output

        result = runner.invoke(cli.get_coauthors, "'Yalchyk'")
        assert result
        assert result.exit_code == 0
        assert result.output == "No coauthors for 'Yalchyk'\n"

    def test_get_common_publications(self):
        runner = CliRunner()

        result = runner.invoke(cli.get_common_publications, "'Birkenbihl Colin' 'Domingo-Fernandez Daniel'")
        assert result
        assert result.exit_code == 0
        assert "33344750" in result.output

        result = runner.invoke(cli.get_common_publications, "'Birkenbihl Colin' 'Yalchyk'")
        assert result
        assert result.exit_code == 0
        assert result.output == "No common publications for 'Birkenbihl Colin' and 'Yalchyk'\n"

    def test_get_common_publications_for_all_coauthors(self):
        runner = CliRunner()

        result = runner.invoke(cli.get_common_publications_for_all_coauthors, "'Birkenbihl Colin'")
        assert result
        assert result.exit_code == 0
        assert "('Domingo-Fernandez Daniel', 'Salimi Yasamin'): ['33344750']" in result.output

        result = runner.invoke(cli.get_common_publications_for_all_coauthors, "'Yalchyk'")
        assert result
        assert result.exit_code == 0
        assert result.output == "No coathors for 'Yalchyk'\n"

    def test_get_common_publications_counts_for_all_coauthors(self):
        runner = CliRunner()

        result = runner.invoke(cli.get_common_publications_counts_for_all_coauthors, "'Birkenbihl Colin'")
        assert result
        assert result.exit_code == 0
        assert "('Domingo-Fernandez Daniel', 'Salimi Yasamin'): 1" in result.output

        result = runner.invoke(cli.get_common_publications_counts_for_all_coauthors, "'Yalchyk'")
        assert result
        assert result.exit_code == 0
        assert result.output == "No coathors for 'Yalchyk'\n"

    def test_generate_graph(self):
        runner = CliRunner()
        TMP_GRAPH_FILE = "graph.png"

        result = runner.invoke(cli.generate_graph, f"'Birkenbihl Colin' graph -f png")
        assert result
        assert result.exit_code == 0
        assert os.path.exists(TMP_GRAPH_FILE)
        os.remove(TMP_GRAPH_FILE)

        result = runner.invoke(cli.generate_graph, "'Yalchyk' graph.png -f png")
        assert result
        assert result.exit_code == 0
        assert result.output == "No coathors for 'Yalchyk'\n"
        assert not os.path.exists(TMP_GRAPH_FILE)
