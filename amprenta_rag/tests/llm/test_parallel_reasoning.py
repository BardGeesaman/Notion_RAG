def test_parallel_query_defaults(monkeypatch):
    # Mock return value of get_available_models to match what parallel_query expects (list of dicts)
    monkeypatch.setattr(pr, "get_available_models", lambda: [{"name": "m1"}, {"name": "m2"}])
    
    # Mock both stages
    monkeypatch.setattr(
        pr, 
        "run_parallel_models", 
        lambda p, models: [{"model": m, "response": "r", "success": True} for m in models]
    )
    monkeypatch.setattr(pr, "synthesize_responses", lambda q, res: "Final")
    
    result = pr.parallel_query("q")
    assert result["question"] == "q"
    assert result["synthesized_answer"] == "Final"
    assert len(result["individual_responses"]) == 2
