{% extends 'markdown.tpl' %}

{% block input %}
<pre><code class="language-python">{{ cell.source }}</code></pre>
{% endblock input %}

{% block stream %}
{% if output.name == "stdout" %}
<pre>{{ output.text }}</pre>
{% endif %}
{% endblock stream %}

{% block data_text scoped %}
<pre>{{ output.data['text/plain'] }}</pre>
{% endblock data_text %}

